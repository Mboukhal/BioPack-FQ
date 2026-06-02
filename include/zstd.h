#pragma once
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

size_t ZSTD_compress(void* dst, size_t dstCapacity, const void* src, size_t srcSize, int compressionLevel);
size_t ZSTD_decompress(void* dst, size_t dstCapacity, const void* src, size_t compressedSize);
unsigned ZSTD_isError(size_t code);
const char* ZSTD_getErrorName(size_t code);
size_t ZSTD_compressBound(size_t srcSize);
unsigned long long ZSTD_getFrameContentSize(const void* src, size_t srcSize);

typedef struct ZSTD_CCtx_s ZSTD_CCtx;
ZSTD_CCtx* ZSTD_createCCtx(void);
size_t ZSTD_freeCCtx(ZSTD_CCtx* cctx);
typedef enum {
    ZSTD_c_nbWorkers=400,
} ZSTD_cParameter;
size_t ZSTD_CCtx_setParameter(ZSTD_CCtx* cctx, int param, int value);
size_t ZSTD_compress2(ZSTD_CCtx* cctx, void* dst, size_t dstCapacity, const void* src, size_t srcSize);

#ifdef __cplusplus
}
#endif
