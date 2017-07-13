using Documenter, MIToS

makedocs(
    format = :html,
    sitename = "MIToS",
    modules = [MIToS],
    pages = [
        "index.md"
    ]
)
