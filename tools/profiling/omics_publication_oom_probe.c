#define _DEFAULT_SOURCE

#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>

int main(void) {
    const size_t allocation_bytes = 512ULL * 1024ULL * 1024ULL;
    const size_t page_bytes = 4096ULL;
    volatile uint8_t *payload = malloc(allocation_bytes);

    if (payload == NULL) {
        return 2;
    }
    for (size_t offset = 0; offset < allocation_bytes; offset += page_bytes) {
        payload[offset] = (uint8_t)(offset / page_bytes);
        if ((offset / page_bytes) % 256ULL == 0ULL) {
            usleep(1000);
        }
    }
    pause();
    return 3;
}
