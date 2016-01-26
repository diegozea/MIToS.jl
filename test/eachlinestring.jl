# using Base.Test
# using MIToS.Utils

print("""

Test eachline for ASCIIString
=============================
""")

let example = ""
  @test_throws ErrorException eachline(example)
end

let example = "Hola"
  @test collect(eachline(example)) == ["Hola"]
end

let example = "Hola\n"
  @test collect(eachline(example)) == ["Hola\n"]
end

let example = "\n"
  @test collect(eachline(example)) == ["\n"]
end

let example = "Hola\nMundo"
  @test collect(eachline(example)) == ["Hola\n", "Mundo"]
end

let example = "Hola\nMundo\n"
  @test collect(eachline(example)) == ["Hola\n", "Mundo\n"]
end

let example = "Hola\nMundo\n\n"
  @test collect(eachline(example)) == ["Hola\n", "Mundo\n", "\n"]
end
