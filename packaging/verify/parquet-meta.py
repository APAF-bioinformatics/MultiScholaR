# run-win-k147 diagnostic (instrumentation): print the Parquet footer of each named file with pyarrow (installed in a
# scratch folder given first). Usage: python parquet-meta.py <pyarrow folder> <file>...
import sys

sys.path.insert(0, sys.argv[1])
import pyarrow
import pyarrow.parquet as pq

print("META pyarrow", pyarrow.__version__)
for path in sys.argv[2:]:
    f = pq.ParquetFile(path)
    m = f.metadata
    print("META ==", path)
    print("META created_by", m.created_by, "| format_version", m.format_version, "| rows", m.num_rows,
          "| row groups", m.num_row_groups, "| columns", m.num_columns, "| serialized_size", m.serialized_size)
    print("META key_value_metadata", m.metadata)
    print("META schema", str(m.schema).replace("\n", " ; "))
    for g in range(m.num_row_groups):
        rg = m.row_group(g)
        print("META row_group", g, "rows", rg.num_rows, "total_byte_size", rg.total_byte_size)
        for c in range(rg.num_columns):
            col = rg.column(c)
            st = col.statistics
            stats = None
            if st is not None:
                stats = dict(has_min_max=st.has_min_max, min=st.min if st.has_min_max else None,
                             max=st.max if st.has_min_max else None, null_count=st.null_count,
                             distinct=st.distinct_count, physical=st.physical_type)
            print("META column", g, c, col.path_in_schema, col.physical_type, "encodings", col.encodings,
                  "compression", col.compression, "compressed", col.total_compressed_size,
                  "uncompressed", col.total_uncompressed_size, "dict_page", col.dictionary_page_offset,
                  "data_page", col.data_page_offset, "stats", stats)
