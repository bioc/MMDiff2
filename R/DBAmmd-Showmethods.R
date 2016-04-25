# Show methods

setMethod("show", "DBAmmd",
          function(object){
            cat(class(object), "instance with",
                numPeaks(object), "peaks and",
                numSamples(object), "samples on\n",
                Genome(object), "genome: \n\n")
            print(Samples(object))
          }
)
