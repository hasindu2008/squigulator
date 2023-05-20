# distutils: language = c
# cython: language_level=3
# cython: profile=True
import sys
import time
import logging
import copy
from libc.stdlib cimport malloc, free
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport strdup
cimport pysigfish
# Import the Python-level symbols of numpy
import numpy as np
# Import the C-level symbols of numpy
cimport numpy as np
# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()


cdef class start:
    # '''
    # Creates a new sigfish object

    # '''
    # cdef sigfish_state_t *state
    # cdef char* REF
    # cdef char* out_paf
    # cdef int NUM_CHANNELS
    # cdef int NUM_THREADS
    # cdef sigfish_read_t *sbatch
    # cdef sigfish_status *status
    # cdef int batch_len
    # cdef object logger
    # cdef char* rid
    # cdef int8_t no_full_ref
    # cdef sigfish_opt_t opt
    # cdef float dtw_cutoff
    # cdef int query_size_sig
    # cdef int query_size_events


    # def __cinit__(self, ref, paf, channels=512, threads=8, dtw_cutoff=70.0, query_size_sig=6000, query_size_events=250, DEBUG=0):
    #     '''
    #     C init
    #     '''
    #     self.state = NULL
    #     self.REF = NULL
    #     self.NUM_CHANNELS = 0
    #     self.NUM_THREADS = 0
    #     self.sbatch = NULL
    #     self.status = NULL
    #     self.batch_len = 0
    #     self.rid = NULL
    #     self.out_paf = NULL
    #     self.no_full_ref = 0
    #     self.dtw_cutoff = 0
    #     self.query_size_sig = 0
    #     self.query_size_events = 0

    #     # sets up logging level/verbosity
    #     self.logger = logging.getLogger(__name__)
    #     if DEBUG == 1:
    #         lev = logging.DEBUG
    #         self.logger.setLevel(logging.DEBUG)
    #     else:
    #         lev = logging.WARNING
    #         self.logger.setLevel(logging.WARNING)

    #     logging.basicConfig(format='%(asctime)s - [%(levelname)s]: %(message)s',
    #                         datefmt='%d-%b-%y %H:%M:%S', level=lev)

    #     print("initiating sigfish", file=sys.stderr)

    #     REF = str.encode(ref)
    #     self.REF = strdup(REF)
    #     PAF = str.encode(paf)
    #     self.out_paf = strdup(PAF)
    #     self.NUM_CHANNELS = channels
    #     self.NUM_THREADS = threads
    #     self.dtw_cutoff = dtw_cutoff
    #     self.query_size_sig = query_size_sig
    #     self.query_size_events = query_size_events

    #     # set up opt struct
    #     self.opt.num_thread = self.NUM_THREADS
    #     self.opt.debug_paf = self.out_paf
    #     self.opt.no_full_ref = self.no_full_ref
    #     self.opt.dtw_cutoff = self.dtw_cutoff
    #     self.opt.query_size_sig = self.query_size_sig
    #     self.opt.query_size_events = self.query_size_events

    #     self.state = init_sigfish(self.REF, self.NUM_CHANNELS, self.opt)
    #     if self.state is NULL:
    #         self.logger.error("Ref '{}' could not be opened and sigfish not initialised".format(ref))

    #     if self.state is NULL:
    #         raise MemoryError()


    # def __init__(self, ref, paf, channels=512, threads=8, dtw_cutoff=70.0, query_size_sig=6000, query_size_events=250, DEBUG=0):
    #     '''
    #     python init
    #     '''
    #     pass

    # def __dealloc__(self):
    #     '''
    #     free memory
    #     '''
    #     if self.out_paf is not NULL:
    #         free(self.out_paf)
    #     # free(self.opt)
    #     # print("free state", file=sys.stderr)
    #     if self.state is not NULL:
    #         free_sigfish(self.state)
    #     # print("free REF", file=sys.stderr)
    #     if self.REF is not NULL:
    #         free(self.REF)



    # def process_batch(self, batch, signal_dtype):
    #     '''
    #     process a batch of of signals
    #     ctypedef struct sigfish_read_t:
    #         int32_t read_number;
    #         int32_t channel;
    #         uint64_t len_raw_signal;
    #         float* raw_signal;

    #     cdef enum sigfish_status:
    #         SIGFISH_MORE = 0,      #more data needed
    #         SIGFISH_REJECT = 1,    #reject the read
    #         SIGFISH_CONT = 2        #continue with the read

    #     readID = read.id
    #     read_number = read.number
    #     chunk_length = read.chunk_length
    #     raw_data = numpy.fromstring(read.raw_data, dtype)
    #     '''
    #     # print("allocating sbatch memory", file=sys.stderr)
    #     status_dic = {}
    #     self.batch_len = len(batch)
    #     if self.batch_len < 1:
    #         self.logger.debug("Batch is of length: {}".format(self.batch_len))
    #         return status_dic

    #     self.sbatch = <sigfish_read_t *> PyMem_Malloc(sizeof(sigfish_read_t)*self.batch_len)
    #     if not self.sbatch:
    #         raise MemoryError()
    #     # print("batch data:", file=sys.stderr)
    #     # for channel, read in batch:
    #     #     print("channel: {}, read_number: {}".format(channel, read.number), file=sys.stderr)
    #     #     break

    #     # print("starting build sbatch for loop", file=sys.stderr)
    #     try:
    #         idx = 0
    #         for channel, read in batch:
    #             self.sbatch[idx].read_number = read.number
    #             rid = str.encode(read.id)
    #             self.sbatch[idx].read_id = strdup(rid)
    #             self.sbatch[idx].channel = channel
    #             self.sbatch[idx].len_raw_signal = read.chunk_length
    #             self.sbatch[idx].raw_signal = <float *> PyMem_Malloc(sizeof(float)*read.chunk_length)
    #             # sig = np.fromstring(read.raw_data, signal_dtype)
    #             # memview = memoryview(sig)
    #             memview = memoryview(np.frombuffer(read.raw_data, signal_dtype))
    #             for i in range(read.chunk_length):
    #                 self.sbatch[idx].raw_signal[i] = memview[i]
    #             idx += 1

    #         # print("calling process_sigfish", file=sys.stderr)
    #         self.status = process_sigfish(self.state, self.sbatch, self.batch_len)
    #         # print("process_sigfish done", file=sys.stderr)

    #         # print("building return", file=sys.stderr)
    #         idx = 0
    #         for channel, read in batch:
    #             status_dic[channel] = (channel, read.number, read.id, self.status[idx], read.raw_data)
    #             idx += 1

    #         # print("returning", file=sys.stderr)
    #         return status_dic
    #     finally:
    #         # print("freeing memory", file=sys.stderr)
    #         # free memory
    #         if self.batch_len > 0:
    #             for i in range(self.batch_len):
    #                 PyMem_Free(self.sbatch[i].raw_signal)
    #                 free(self.sbatch[i].read_id)
    #                 # for j in range(self.sbatch[i].len_raw_signal):
    #                     # free(self.sbatch[i].raw_signal[j])
    #         # print("free sbatch", file=sys.stderr)
    #         PyMem_Free(self.sbatch)
    #         # print("free status", file=sys.stderr)
    #         free(self.status)
    #         # print("free rid", file=sys.stderr)
    #         # free(self.rid)
    #         # free(self.rid)


