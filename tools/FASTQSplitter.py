import sys
import zlib
import subprocess
import os
import uuid
import shutil
import io
import logging
import traceback

logger = logging.getLogger('chunkwise_logger')


# constant
##############################################

SNZIP = '/usr/local/bin/snzip'
HADOOP = '/usr/local/hadoop/bin/hadoop'
TMP_DIR = '/tmp'
BUFFER_READ_SIZE = 16384
GZ_SPLIT_CHUNK_SIZE = 256 * 1024 * 1024
SUBPROCESS_NUM = 3
PLAIN_TEXT_LENGTH_LIMIT = 600 * 1024 * 1024
HDFS_NAME_TEMPLATE = '{}/chunk_{}.fastq.snappy'
FS_NAME_TEMPLATE = '{}/tmp_chunk-{}.fq'
FS_NAME_PAIRED_TEMPLATE = '{}/tmp_chunk-{}-{}.fq'
CONTENT_BUFFER_SIZE = 10240
SNZIP_N_UPLOAD_RETRY = 20
PAIRED_TWO_SEARCH_SCOPE = 1 * 1024 * 1024
PAIR_TWO_SEARCH_BUF_SIZE = 10 * 1024

# class
##############################################


class DatasetsUploadError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


class GzChunkFileWrapper:
    def __init__(self, source, is_dir=True, input_chunk_size=GZ_SPLIT_CHUNK_SIZE, handle_buff_size=1024 * 1024):
        self.decompressor = zlib.decompressobj(zlib.MAX_WBITS | 16)
        self.input_chunk_size = input_chunk_size
        self.index = 0
        self.source = source
        self.is_dir = is_dir

        if not os.path.isdir(self.source) and is_dir:
            logger.error("provided source not a directory {} when is_dir==True".format(self.source))
            raise Exception("provided source not a directory {} when is_dir==True".format(self.source))
        if os.path.isdir(self.source) and not is_dir:
            logger.error("provided source not a file {} when is_dir!=True".format(self.source))
            raise Exception("provided source not a file {} when is_dir!=True".format(self.source))

        self.path_list, self.path_list_size, self.list_last_index = self.init_input_file_list()
        self.source_handle, self.buff, self.handle = self.init_file_handle(self.input_chunk_size, handle_buff_size)
        logger.info('path_list {}; path_list_length {}'.format(self.path_list, self.path_list_size))

    def __del__(self):
        self.buff.close()
        self.handle.close()

    def init_input_file_list(self):
        if self.is_dir:
            tmp = [f for f in os.listdir(self.source)]
            logger.debug('input chunk list -- {}'.format(tmp))
            tmp.sort(key=int)
            chunk_list = ["{}/{}".format(self.source, f) for f in tmp]
        else:
            chunk_list = [self.source]
        return chunk_list, len(chunk_list), len(chunk_list) - 1

    def init_file_handle(self, chunk_size, handle_buff_size):
        f = open(self.path_list[self.index], "rb")
        buff = io.BytesIO(f.read(chunk_size))
        return f, buff, io.BufferedReader(buff, handle_buff_size)

    def read(self, length=BUFFER_READ_SIZE):
        try:
            total = 0
            buffer = self.handle.read(length)
            # return length < read length => self.buff is consumed, and need to load next chunk to self.buff
            if len(buffer) < length:
                logger.debug('buffer.len -- {}'.format(len(buffer)))
                self.append_next_chunk()
                buffer += self.handle.read(length)

            # decompress buffer (gz.compressed) to outstr (binary string)
            outstr = self.decompressor.decompress(buffer)
            total += len(outstr)

            # loop through all the unused_data, and append decompressed portion into outstr
            while self.decompressor.unused_data != b'':
                unused_data = self.decompressor.unused_data
                self.decompressor = zlib.decompressobj(zlib.MAX_WBITS | 16)
                tmp = self.decompressor.decompress(unused_data)
                total += len(tmp)
                outstr += tmp
        except:
            logger.error('input decompress error -- {}'.format(traceback.format_exc()))
            raise DatasetsUploadError('input decompress error -- {}'.format(traceback.format_exc()))

        return len(outstr) == 0, outstr

    def append_next_chunk(self):
        # check whether source_handle has been consumed
        if self.index > self.list_last_index:
            return

        content = self.source_handle.read(self.input_chunk_size)

        if not content:
            # remove consumed chunk file
            self.source_handle.close()
            if self.is_dir:
                os.unlink(self.path_list[self.index])

            # about to append next chunk
            self.index += 1

            if self.index > self.list_last_index:
                return

            self.source_handle = open(self.path_list[self.index], 'rb')

            # append chunk routine
            remain_buf = self.buff.read(sys.getsizeof(self.buff) - self.buff.tell())
            self.buff.truncate(0)
            self.buff.seek(0)
            self.buff.write(remain_buf)
            self.buff.write(self.source_handle.read(self.input_chunk_size))
            logger.debug("__appended chunk No. {}, with buf length -- {}".format(self.index, sys.getsizeof(self.buff)))
            self.buff.seek(0)
        else:
            remain_buf = self.buff.read(sys.getsizeof(self.buff) - self.buff.tell())
            self.buff.truncate(0)
            self.buff.seek(0)
            self.buff.write(remain_buf)
            self.buff.write(content)
            logger.debug("appending next source file content, with buf length -- {}".format(sys.getsizeof(self.buff)))
            self.buff.seek(0)


class FastqUploader:
    def __init__(self, file_dir_list, dest_path,
                 fastq_txt_chunk_size=PLAIN_TEXT_LENGTH_LIMIT,
                 input_chunk_size=GZ_SPLIT_CHUNK_SIZE):
        self.upload_retried_count = 0
        self.source_handle = []
        logger.info(file_dir_list)
        logger.info(dest_path)
        for item in file_dir_list:
            if os.path.isdir(item):
                self.source_handle.append(GzChunkFileWrapper(item, True, input_chunk_size))
            elif os.path.isfile(item):
                self.source_handle.append(GzChunkFileWrapper(item, False, input_chunk_size))
            else:
                logger.error('input source {} is not a directory nor a file'.format(item))
                raise DatasetsUploadError('input source {} is not a directory nor a file'.format(item))

        self.result_dest = dest_path
        self.tmp_dir = '{}/{}'.format(TMP_DIR, uuid.uuid4())
        try:
            os.mkdir(self.tmp_dir)
        except FileExistsError:
            logger.error("{} -- already exist".format(self.tmp_dir))
        self.init_dest_dir()
        # for single-end file, self.fastq_txt_chunk_size = fastq_txt_chunk_size
        self.fastq_txt_chunk_size = int(fastq_txt_chunk_size / len(self.source_handle))
        logger.debug('source dir: {}, dest dir: {}, chunk fastq plain text size {}'.format(
            file_dir_list, format(dest_path), self.fastq_txt_chunk_size))

    def init_dest_dir(self):
        return_code = subprocess.call([HADOOP, 'fs', '-test', '-e', self.result_dest])
        if return_code != 0:
            subprocess.check_call([HADOOP, 'fs', '-mkdir', self.result_dest])
        logger.info('result dest dir -- {}'.format(self.result_dest))

    @staticmethod
    def find_split_point(shortfall, content):
        content_size = len(content)
        idx1 = shortfall
        while idx1 < content_size:
            if content[idx1] == ord('\n') and content[idx1 + 1] == ord('@'):
                idx2 = idx1 + 1
                while idx2 < content_size:
                    if content[idx2] == ord('\n'):
                        if content[idx2 + 1] == ord('@'):
                            return idx2 + 1
                        else:
                            return idx1 + 1
                    idx2 += 1
            idx1 += 1
        logger.error("cannot find proper position for partitioning fastq files")
        return -1

    def process_file(self, chunk_file_path, chunk_index, pid):
        # compress and upload file
        cmd = '{0} -k -t hadoop-snappy {1} && {4} fs -put {2} {3}'.format(
            SNZIP,
            chunk_file_path,
            '{}.snappy'.format(chunk_file_path,),
            HDFS_NAME_TEMPLATE.format(self.result_dest, str(chunk_index).zfill(5)),
            HADOOP)
        logger.debug(cmd)
        pid.append((subprocess.Popen(cmd, shell=True), cmd, [chunk_file_path, '{}.snappy'.format(chunk_file_path)]))
        if len(pid) > SUBPROCESS_NUM:
            return self.join_subprocess(pid)
        return pid

    def join_subprocess(self, pid):
        retry_pid = []
        for item in pid:
            return_code = item[0].wait()
            if return_code != 0:
                if self.upload_retried_count < SNZIP_N_UPLOAD_RETRY:
                    self.upload_retried_count += 1
                    logger.error('retry {} with command {}'.format(self.upload_retried_count, item[1]))
                    retry_pid.append((subprocess.Popen(item[1], shell=True), item[1], item[2]))
                else:
                    logger.error('compress process error failed {} time'.format(SNZIP_N_UPLOAD_RETRY))
                    raise DatasetsUploadError('compress process error failed {} time'.format(SNZIP_N_UPLOAD_RETRY))
            else:
                for f in item[2]:
                    try:
                        os.remove(f)
                    except:
                        logger.error('failed to remove {}'.format(f))
        pid = retry_pid
        return pid

    def run(self):
        try:
            chunk_index = 0
            written_length = 0
            chunk_file_path = FS_NAME_TEMPLATE.format(self.tmp_dir, chunk_index)
            chunk_f = open(chunk_file_path, 'wb')
            pid = []
            while True:
                eof, content = self.source_handle[-1].read()
                # shortfall_length = self.fastq_txt_chunk_size - written_length
                # len(content) > shortfall_length => reach current chunk length threshold, and find the
                # chunk-segmentation at the start of the 1st read right after content[:shortfall_length]
                shortfall_length = self.fastq_txt_chunk_size - written_length

                # CONTENT_BUFFER_SIZE is a buffer size of content, which is set to be 10240,
                # to guarantee that content is at least 10240 characters longer than shortfall_length for the following
                # split point search operation.
                # CONTENT_BUFFER_SIZE should be configurable based on size of reads.
                # For Illumina short reads scenario, CONTENT_BUFFER_SIZE can be set to 512 since size of Illumina
                # short reads rarely exceeds 300 byte.
                if len(content) > shortfall_length + CONTENT_BUFFER_SIZE:
                    if shortfall_length < 0:
                        shortfall_length = 0
                    split_point = self.find_split_point(shortfall_length, content)
                    written_length += split_point
                    chunk_f.write(content[:split_point])
                    logger.debug('done {} writing, file length {}'.format(chunk_file_path, written_length))
                    chunk_f.close()
                    pid = self.process_file(chunk_file_path, chunk_index, pid)

                    chunk_index += 1
                    chunk_file_path = FS_NAME_TEMPLATE.format(self.tmp_dir, chunk_index)
                    chunk_f = open(chunk_file_path, 'wb')
                    chunk_f.write(content[split_point:])
                    written_length = len(content) - split_point

                else:
                    written_length += len(content)
                    chunk_f.write(content)

                if eof:
                    # eof chunk writing
                    written_length += len(content)
                    chunk_f.write(content)
                    logger.debug('done eof chunk - {} writing, file length {}'.format(chunk_file_path, written_length))
                    chunk_f.close()
                    self.process_file(chunk_file_path, chunk_index, pid)
                    self.join_subprocess(pid)
                    break

        finally:
            pass
            # clean tmp files and source chunks
            shutil.rmtree(self.tmp_dir, ignore_errors=True)


def main(args):
    uploader = FastqUploader([sys.argv[1]], sys.argv[2])
    uploader.run()


if __name__ == '__main__':
    main(sys.argv)

