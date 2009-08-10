#########################################################################
#                                                                       #
# GenK: Genomic Ripley's K Using the Network Method                     #
#                                                                       #
# Author: Jamie Estill                                                  #
# Contact: jestill@uga.edu                                              #
# Started: 10/27/2004                                                   #
# Last Updated: 9/21/2005                                               #
#                                                                       #
# Usage:                                                                #
# GenKh(h, L, P)                                                        #
#       h = CriticalDistance
#           This is the lag distance that will be used.                 #
#       L = Length of the entire genomic entire being considered        #
#       P = Vector of point location of features being tested           #
# Exampe Usage:                                                         #
# NetKh(100000, 3000000000, ArabGenes)                                  #
#                                                                       #
# This function will implement a Network Ripley's K as described in:    #
# Yamada and Thill. 2004. Journal of Transport Geography. 12: 149-158   #
# which is a slight modification of                                     #
# Okabe and Yamada. 2001. Geographic Analysis. 33(3):271-290            #
# This genomic network K will take as its input the set of positions    #
# of the features that are of interest.
#                                                                       #
#########################################################################

#########################################################################
#                                                                       #
# Function to calculate K for a given critical distance h               #
# This is in the simplified summation approach.                         #
# The needs to be converted to a matrix method.                         #
#                                                                       #
#########################################################################

GenKh <- function(CritDistance, Length, Features)

{ # Begin of GenKh function

  #print (paste("Critical Distance:",CritDistance))
  #print (paste("Length:",Length))
  #print (paste("Features:",Features))
  #print (paste("Feature01:",Features[1]))

  # This gave a out of bounds error
  #print (paste("Feature02:",Features[1,2]))


  NumFeatures = length(Features)          # Calculate the number of features in the list

  Density = Length/(NumFeatures*(NumFeatures-1))
  Sum = 0                                    # Set the sum to zero for starting

  for (i in 1:NumFeatures)
  {   # Begin of i
      for (j in 1:NumFeatures)
      { # Begin of j
        if (i != j) # If i is not equal to j ... in other words do not compare values to themselves
        { # Begin of if i not equal j

          # Following print statements are for debug purposes
          #print (paste("I Value:",i,"J Value",j))
          #print (paste("Compare",Features[i],Features[j]))
          #print (paste("Sum",Sum))
          EvalDistance = abs(Features[i]-Features[j])
          # The following did not work so I changed it to the above
          # EvalDistance = abs(Features[i,1]-Features[j,1])
          if (EvalDistance <= CritDistance)
          {
           #print (paste("DING","Compare",Features[i],Features[j]))
           Sum = Sum + 1
          }
        } # End of i not equl to j
      } # End of j
  }   # End of for i
  
print (paste("Sum",Sum))
print (paste("Density",Density))
Kh = Density*Sum    # calcluate Kh
Kh                  # this will return the calculate Kh value
} # End of GenKh function



WorkKh <- function(CritDistance, Length, Features)

{ # Begin of GenKh function

  #print (c("Critical Distance",CritDistance))
  #print (c("Length",Length))
  #print (c("Features",Features))

  NumFeatures = length(Features)          # Calculate the number of features in the list
  #Density = Length/(NumFeatures*(NumFeatures-1))
  #Sum = 0                                    # Set the sum to zero for starting
  #
  #for (i in 1:NumFeatures)
  #{   # Begin of i
  #    for (j in 1:NumFeatures)
  #    { # Begin of j
  #      if (i != j) # If i is not equal to j ... in other words do not compare values to themselves
  #      { # Begin of if i not equal j
  #        EvalDistance = abs(Features[i,1]-Features[j,1])
  #        if (EvalDistance <= CritDistance)
  #        {
  #         Sum = Sum + 1
  #        }
  #      } # End of i not equl to j
  #    } # End of j
  #}   # End of for i

#Kh = Density*Sum    # calcluate Kh
#Kh                  # this will return the calculate Kh value

NumFeatures
} # End of GenKh function



