setGeneric("foo", function(x) standardGeneric("foo"))
setGeneric("foo", function(y) standardGeneric("foo"))
setClass("A")
setMethod("foo", "A", function(x) "A")
