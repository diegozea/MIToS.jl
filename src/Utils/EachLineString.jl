"Create an iterable object that will yield each line from a stream **or string**."
lineiterator(string::String) = eachline(IOBuffer(string))
lineiterator(stream::IO)     = eachline(stream)

using Requests

# allow_redirects = false
# headers = Dict("User-Agent" => "Mozilla/5.0 (compatible; MSIE 7.01; Windows NT 5.0)")
function download_file(url::AbstractString, filename::AbstractString; kargs...)
    stream = Requests.get_streaming(url; kargs...)
    open(filename, "w") do fh
        while !eof(stream)
            write(fh, readavailable(stream))
        end
    end
    filename
end

download_file(url::AbstractString; kargs...) = download_file(url, tempname(); kargs...)
