#include <stdlib.h>

void* (*stbiMallocPtr)(size_t size) = NULL;
void* (*stbiReallocPtr)(void* ptr, size_t size) = NULL;
void (*stbiFreePtr)(void* ptr) = NULL;

#define STBI_MALLOC(size) stbiMallocPtr(size)
#define STBI_REALLOC(ptr, size) stbiReallocPtr(ptr, size)
#define STBI_FREE(ptr) stbiFreePtr(ptr)

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

void* (*stbirMallocPtr)(size_t size, void* context) = NULL;
void (*stbirFreePtr)(void* ptr, void* context) = NULL;

#define STBIR_MALLOC(size, context) stbirMallocPtr(size, context)
#define STBIR_FREE(ptr, context) stbirFreePtr(ptr, context)

#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "stb_image_resize2.h"

void* (*stbiwMallocPtr)(size_t size) = NULL;
void* (*stbiwReallocPtr)(void* ptr, size_t size) = NULL;
void (*stbiwFreePtr)(void* ptr) = NULL;

#define STBIW_MALLOC(size) stbiwMallocPtr(size)
#define STBIW_REALLOC(ptr, size) stbiwReallocPtr(ptr, size)
#define STBIW_FREE(ptr) stbiwFreePtr(ptr)

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
