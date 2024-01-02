#OK I'm jut going to have to figure out the parital R-squared myself.
#start with a glmmTMP obje




PartialR2 = function(mod){
  #start with the count model
  terms = attr(mod$modelInfo$terms$cond$fixed, "term.labels")
  termszi = attr(mod$modelInfo$terms$zi$fixed, "term.labels")
  
  #set up data frame to hold the results
  out = data.frame(terms = c(terms, termszi), 
                   paramtype = c(rep("count", length(terms)), rep("zi",length(termszi))),
                   R2 = NA)
  Form = mod$modelInfo$terms$cond$fixed
  FormZI = mod$modelInfo$terms$zi$fixed
  for(i in 1:length(terms)) {
    #set up reduced model
    term = noquote(terms[i])
    FormX = update(Form, ~ . + eval(term))
    reduced = glmmTMB(FormX, 
                      zi = ~ termszi + (1|Region) + (1|Year),data = mod$frame,
                      family= "nbinom2")
    perf = performance(reduced)
out$
    
  }
}

X = model.matrix(mod)
drop = which(colnames(X) == 'education:typeprof')
X1 = X[,-1]
newfit = lm(presitge ~ X1 - 1)

library(glmmTMB)
library(performance)
fm <- glmmTMB(y ~ Base+trt + Age + Visit + (1|subject),
              data=epil2, family=nbinom2)

PartialR2 = function(mod){
  #Pull the terms out of the model
  terms = attr(mod$modelInfo$terms$cond$fixed, "term.labels")
  #set up data frame to hold the results
  out = data.frame(Terms = c(terms),
                   R2 = NA)
  #Pull the formula out of the model
  Form = mod$modelInfo$terms$cond$fixed
  for(i in 1:length(terms)) {
    #set up reduced model. THIS IS THE BIT THAT ISN"T WORKING!!!
    FormX = update(Form, ~ . - terms[i])
    reduced = glmmTMB(FormX, data = mod$frame,
                      family= "nbinom2")
    perf = performance(reduced)
      out$Terms[i] = terms[i]
      out$R2 = perf$R2_marginal
  }
  return(out)
}


###############################################
#try a different way
#start with a dredge object

Rs = function(d){
  foo = dplyr::select(as.data.frame(d), `cond((Int))`:`df`)
  foo$nterms = apply(!is.na(foo), 1, sum)-3
test = dplyr::select(foo, -Rs.Rsqc, -Rs.Rsqm, -df, -nterms)
Nterms = foo$nterms[1]
terms = names(test[1,which(!is.na(test[1,]))])

R2best = foo$Rs.Rsqc[1]
foosub = dplyr::filter(foo, nterms == (Nterms -1)) %>%
  mutate(rpart = R2best - Rs.Rsqc,)

for(i in 1:nrow(foosub)) {
  foosub2 = foosub[i,]
  term = names(foosub2[,which(is.na(foosub2))])
  term %in% terms
  term = terms[which(names(foosub[i]) %in% terms)]
  foosub$term[i] = term 
}

}

#I GIVE UP
