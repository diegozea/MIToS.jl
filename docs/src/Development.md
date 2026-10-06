# Development

This section describes how to extend the MIToS API.

## [MSA](@id Development-MSA)

### [Extending eachmsa](@id Extending-eachmsa)

To make [`eachmsa`](@ref MIToS.MSA.eachmsa) work with another `MSAFormat`:

  - Define `MIToS.MSA.support_eachmsa(::Type{MyFormat}) = true` to declare support.
  - Implement [`hasnextmsa`](@ref MIToS.MSA.hasnextmsa) to check for another alignment,
    leaving its beginning available to read.
  - Implement [`parse_file`](@ref MIToS.Utils.parse_file) to read one alignment at a time,
    accepting the requested output type and reading options.

Both reading methods receive the same buffered input and should leave it open.
See [`hasnextmsa`](@ref MIToS.MSA.hasnextmsa) for using a header pattern, handling
unexpected input, and putting a line back for the next alignment.
[`read_file`](@ref MIToS.Utils.read_file) also uses this interface to warn when more
alignments follow the first one.
