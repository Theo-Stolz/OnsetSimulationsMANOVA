permtfdist <- function(cond1, cond2, cond1.grad, cond2.grad, Nf, Nt, nboot = 2000){
  
  st2 <- Nt+1 # index to start condition 2 
  en2 <- Nt+Nt # index to end condition 2
  all.amp <- Rfast::transpose(cbind(cond1, cond2)) # combine trials, then transpose -> trials x time points
  all.grad <- cbind(cond1.grad, cond2.grad) # combine trials only
  
  perm.t.amp <- matrix(data = 0, nrow = nboot, ncol = Nf)
  perm.f <- matrix(data = 0, nrow = nboot, ncol = Nf)
  
  for(B in 1:nboot){
    
    # permutation indices
    perm.ind <- sample(1:en2, size = en2, replace = FALSE)
    
    # t-test on amplitudes
    perm.amp <- all.amp[perm.ind,]
    perm.amp1 <- perm.amp[1:Nt,]
    perm.amp2 <- perm.amp[st2:en2,]
    perm.t.amp[B,] <- Rfast::ttests(perm.amp1, perm.amp2, paired=FALSE)[,1]    
    
    # MANOVA on amplitudes and gradients
    perm.grad <- all.grad[,perm.ind]
    perm.grad1 <- perm.grad[,1:Nt]
    perm.grad2 <- perm.grad[,st2:en2]
    for(F in 2:Nf){
      perm.f[B,F] <- hotelling.stat(cbind(perm.amp1[,F],perm.grad1[F,]), cbind(perm.amp2[,F],perm.grad2[F,]),
                                    shrinkage=FALSE, var.equal=FALSE)$statistic
    }
  }
  list(t.amp = perm.t.amp, f.amp.grad = perm.f)
}