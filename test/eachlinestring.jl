#using Base.Test
#using MIToS.Utils

print("""

Test eachline for ASCIIString
=============================
""")

let example = "Hola"
  @test collect(lineiterator(example)) == ["Hola"]
end

let example = "Hola\n"
  @test collect(lineiterator(example)) == ["Hola\n"]
end

let example = "\n"
  @test collect(lineiterator(example)) == ["\n"]
end

let example = "Hola\nMundo"
  @test collect(lineiterator(example)) == ["Hola\n", "Mundo"]
end

let example = "Hola\nMundo\n"
  @test collect(lineiterator(example)) == ["Hola\n", "Mundo\n"]
end

let example = "Hola\nMundo\n\n"
  @test collect(lineiterator(example)) == ["Hola\n", "Mundo\n", "\n"]
end
