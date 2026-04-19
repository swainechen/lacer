## 2025-05-14 - Buffer Overflows and Memory Leaks in lacepr
**Vulnerability:** Multiple buffer overflows via `sscanf`, `strcpy`, and `strcat`, plus a major memory leak in the BAM processing loop.
**Learning:** The C component `lacepr` used unsafe string functions without length checks when parsing recalibration files and BAM tags. It also allocated memory for strings on every read/record without freeing them.
**Prevention:** Always use width limits in `sscanf` (e.g., `%255s`), use `strncpy` with explicit null-termination, and avoid redundant heap allocations in high-frequency loops by pointing to existing buffers or reusing memory.

## 2025-05-15 - Command Injection via two-argument open() in Perl
**Vulnerability:** Use of two-argument `open()` with user-supplied filenames in `lacer.pl`.
**Learning:** In Perl, the two-argument `open(FH, $filename)` can execute commands if `$filename` starts with a pipe (`|`). This is a classic command injection vector.
**Prevention:** Always use the three-argument form `open(my $fh, '<', $filename)` and lexical filehandles. This explicitly separates the open mode from the filename, preventing injection.

## 2025-05-16 - NULL Pointer Dereference and Insecure realloc in lacepr
**Vulnerability:** Use of uninitialized pointers in `access()` and failure to update caller's pointer after `realloc()` in `read_recal`.
**Learning:** Command-line argument pointers in C were passed to `access()` without initialization or NULL checks, leading to crashes if arguments were missing. In `read_recal`, `realloc` was used on a pointer passed by value, meaning the caller still held the potentially freed original pointer.
**Prevention:** Initialize all pointers to `NULL`, validate them before use, and pass pointers-to-pointers (e.g., `recal_t **data_ptr`) to functions that resize memory to ensure the caller's reference is updated.

## 2025-05-17 - OOB Access in Cycle Indexing and Missing Memory Safety in lacepr
**Vulnerability:** Heap buffer overflow when accessing `Cycle` and `cache` arrays with indices from `get_cycle_index`, plus multiple missing NULL checks for `malloc` and `samopen`.
**Learning:** The `Cycle` array was sized to `MAX_CYCLE + 1`, but `get_cycle_index` could return values up to `2 * MAX_CYCLE`. This mismatch caused out-of-bounds writes. Also, `cache` was indexed by quality up to `MAX_Q` but only sized for `MAX_REASONABLE_Q`.
**Prevention:** Ensure array dimensions fully cover the range of all potential index transformation functions. Implement defensive bounds checks at the point of access in high-frequency functions and strictly validate all memory allocations and library resource handles (like `HTSlib`'s `samopen`).

## 2025-05-18 - Sticky Variable State and Missing Library Call Validation in lacepr
**Vulnerability:** Data integrity risk from "sticky" read group indices in BAM processing and potential NULL pointer dereference if `gzopen` fails.
**Learning:** In high-frequency record processing loops, variables like `rg_index` that are conditionally updated can retain state from previous records if not explicitly reset. Additionally, external library calls for file I/O (like `htslib`'s `gzopen`) were assumed to always succeed.
**Prevention:** Always reset per-record state at the beginning of processing loops. Strictly validate the return values of all library-provided resource acquisition functions (e.g., `gzopen`, `samopen`, `malloc`) before usage.

## 2025-05-19 - OOB Array Access via Malformed Recalibration Files in lacepr
**Vulnerability:** Out-of-bounds (OOB) memory access when parsing GATK-style recalibration tables in `read_recal`, where `rgindex`, `qual`, and `covindex` were used as array indices without bounds checking.
**Learning:** External data from files (even those typically produced by trusted tools like `lacer.pl`) must be treated as untrusted. Indices derived from `sscanf` can exceed array dimensions, leading to heap/stack corruption.
**Prevention:** Always implement explicit bounds checks (e.g., `if (index >= 0 && index < MAX_SIZE)`) immediately after parsing and before using the index to access any array or allocated memory.
