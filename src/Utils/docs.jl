# ## Documentation helpers
#
# This functions are used to generate IJulia documentation notebooks for MIToS.
#
# MIToS docstrings don't need to have titles,
# because the titles are displayed by this macros:
#
# User defines the main titles (h1)
# 	- Modules (h2)
#		- Data types (h3)
#	 		- Functions/Methods (h4)
#			- Constants (h4)

"This macro displays the documentation of a module"
macro moduledoc(m)
    quote
        Base.display( Base.Markdown.@md_str($(Base.string("## ", m))) )     # h2
        Base.display( Base.Docs.@repl( $(esc(m)) ) )
        nothing
    end
end

"This macro displays the documentation of a data type"
macro typedoc(t)
    name = Base.split(Base.string(t), '.')[end]
    quote
        Base.display( Base.Markdown.@md_str($(Base.string("### ", name))) ) # h3
        Base.display( Base.Docs.typesummary( $(esc(t)) ) )
        Base.display( Base.Docs.@repl( $(esc(t)) ) )
        nothing
    end
end

"This macro displays the documentation of a function"
macro functiondoc(f)
    name = Base.split(Base.string(f), '.')[end]
    quote
        Base.display( Base.Markdown.@md_str($(Base.string("#### ", name))) )# h4
        Base.display( Base.Docs.@repl( $(esc(f)) ) )
        nothing
    end
end

"This macro displays the documentation of a constant"
macro constantdoc(c)
    name = Base.split(Base.string(c), '.')[end]
    quote
        Base.display( Base.Markdown.@md_str($(Base.string("#### ", name))) )# h4
        Base.display( $(esc(c)) )
        Base.display( Base.Docs.@repl( $(esc(c)) ) )
        nothing
    end
end
