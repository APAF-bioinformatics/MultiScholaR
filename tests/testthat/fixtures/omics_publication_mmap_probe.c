#define _GNU_SOURCE

#include <errno.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

static int write_text(const char *path, const char *text) {
    FILE *stream = fopen(path, "w");
    if (stream == NULL) {
        return 1;
    }
    if (fputs(text, stream) == EOF || fclose(stream) != 0) {
        return 1;
    }
    return 0;
}

static double monotonic_seconds(void) {
    struct timespec value;
    if (clock_gettime(CLOCK_BOOTTIME, &value) != 0) {
        return -1.0;
    }
    return (double)value.tv_sec + (double)value.tv_nsec / 1000000000.0;
}

static int write_resources(const char *run_dir) {
    char path[4096];
    const char *state =
        "{\n"
        "  \"arrow_pool_bytes\": 0,\n"
        "  \"duckdb_memory_bytes\": 0,\n"
        "  \"duckdb_spill_bytes\": 0,\n"
        "  \"duckdb_connections\": 0,\n"
        "  \"duckdb_results\": 0,\n"
        "  \"duckdb_prepared_statements\": 0,\n"
        "  \"temporary_paths\": 0,\n"
        "  \"cache_entries\": 0,\n"
        "  \"active_tasks\": 0,\n"
        "  \"observers\": 0,\n"
        "  \"native_resources\": 0\n"
        "}";
    char json[16384];
    snprintf(
        json,
        sizeof(json),
        "{\n"
        "  \"schema\": \"multischolar.omics_publication_worker_resources\",\n"
        "  \"schema_version\": \"1.0.0\",\n"
        "  \"high_water\": %s,\n"
        "  \"retained\": %s,\n"
        "  \"terminal\": %s\n"
        "}\n",
        state,
        state,
        state
    );
    snprintf(path, sizeof(path), "%s/worker-resources.json", run_dir);
    return write_text(path, json);
}

static int write_retention_state(const char *run_dir) {
    char path[4096];
    char json[8192];
    double settled = monotonic_seconds();
    if (settled < 0.0) {
        return 1;
    }
    snprintf(
        json,
        sizeof(json),
        "{\n"
        "  \"active_tasks\": 0,\n"
        "  \"open_queries\": 0,\n"
        "  \"open_writers\": 0,\n"
        "  \"open_leases\": 0,\n"
        "  \"open_connections\": 0,\n"
        "  \"open_results\": 0,\n"
        "  \"active_child_processes\": 0,\n"
        "  \"arrow_pool_bytes\": 0,\n"
        "  \"duckdb_memory_bytes\": 0,\n"
        "  \"duckdb_spill_bytes\": 0,\n"
        "  \"duckdb_prepared_statements\": 0,\n"
        "  \"temporary_paths\": 0,\n"
        "  \"cache_entries\": 0,\n"
        "  \"observers\": 0,\n"
        "  \"native_resources\": 0,\n"
        "  \"retained_dwell_seconds\": 5,\n"
        "  \"retention_acknowledgement\": \"fifo_v1\",\n"
        "  \"settled_monotonic_seconds\": %.9f\n"
        "}\n",
        settled
    );
    snprintf(path, sizeof(path), "%s/retention-state.json", run_dir);
    return write_text(path, json);
}

int main(int argc, char **argv) {
    const size_t mapping_bytes = 32ULL * 1024ULL * 1024ULL;
    const size_t page_bytes = 4096ULL;
    char data_path[4096];
    char fifo_path[4096];
    char marker_path[4096];
    int data_fd;
    int fifo_fd;
    volatile uint8_t *mapping;
    uint8_t token;

    if (argc != 2) {
        return 2;
    }
    snprintf(data_path, sizeof(data_path), "%s/mmap-page-cache.bin", argv[1]);
    snprintf(fifo_path, sizeof(fifo_path), "%s/retention-sampled.fifo", argv[1]);
    snprintf(marker_path, sizeof(marker_path), "%s/retention-ready", argv[1]);

    data_fd = open(data_path, O_CREAT | O_TRUNC | O_RDWR, 0600);
    if (data_fd < 0 || ftruncate(data_fd, (off_t)mapping_bytes) != 0) {
        return 3;
    }
    mapping = mmap(NULL, mapping_bytes, PROT_READ | PROT_WRITE, MAP_SHARED, data_fd, 0);
    if (mapping == MAP_FAILED) {
        return 4;
    }
    for (size_t offset = 0; offset < mapping_bytes; offset += page_bytes) {
        mapping[offset] = (uint8_t)(offset / page_bytes);
    }
    if (msync((void *)mapping, mapping_bytes, MS_SYNC) != 0) {
        return 5;
    }
    if (mkfifo(fifo_path, 0600) != 0 || write_resources(argv[1]) != 0 ||
        write_retention_state(argv[1]) != 0 ||
        write_text(marker_path, "ready\n") != 0) {
        return 6;
    }

    fifo_fd = open(fifo_path, O_RDONLY);
    if (fifo_fd < 0 || read(fifo_fd, &token, 1) != 1 || token != 1) {
        return 7;
    }
    close(fifo_fd);
    munmap((void *)mapping, mapping_bytes);
    close(data_fd);
    return 0;
}
