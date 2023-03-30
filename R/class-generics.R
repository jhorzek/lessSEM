# we create new generics so that we won't have to import stats which causes
# a namespace clash with lavaan

setGeneric("coef", function(object, ...){
  standardGeneric("coef")
})

setGeneric("plot", function(x, ...){
  standardGeneric("plot")
})