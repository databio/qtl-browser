"""One adapter per data source (pipeline/CONTRACT.md).

An adapter knows one source's file names, column spellings and allele convention, and turns them
into the five contract tables under `data/derived/_tables/<experiment_id>/`. Nothing below that
directory knows where the numbers came from, and nothing above it knows the qtlb format.
"""
