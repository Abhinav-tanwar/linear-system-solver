rref <- function(A,tol=1e-10) {
 A <- as.matrix(A); row <- 1L
 for(col in seq_len(ncol(A))) {
  if(row>nrow(A)) break
  pivot <- row-1L+which.max(abs(A[row:nrow(A),col]))
  if(abs(A[pivot,col])<=tol) next
  A[c(row,pivot),] <- A[c(pivot,row),]; A[row,] <- A[row,]/A[row,col]
  for(i in seq_len(nrow(A))) if(i!=row) A[i,] <- A[i,]-A[i,col]*A[row,]
  row <- row+1L
 }
 A[abs(A)<tol] <- 0; A
}
run_analysis <- function(root,out) {
 set.seed(530301); A <- matrix(sample(-10:10,100,replace=TRUE),10); b <- runif(10,-10,10)
 if(qr(A)$rank<ncol(A)) stop("Generated matrix is singular")
 # Inspect matrix structure without changing the coefficient matrix.
 structure <- list(first_rows=A[1:3,,drop=FALSE],element_2_3=A[2,3],diagonal=diag(A),upper=A*upper.tri(A,diag=TRUE),lower=A*lower.tri(A,diag=TRUE),augmented=cbind(A,b))
 for(name in names(structure)) write.csv(structure[[name]],file.path(out,paste0(name,".csv")),row.names=FALSE)
 B<-matrix(sample(11:20,100,replace=TRUE),10)
 identity_errors<-c(sum=max(abs(t(A+B)-(t(A)+t(B)))),product=max(abs(t(A%*%B)-t(B)%*%t(A))),double=max(abs(t(t(A))-A)))
 set.seed(530311); normal_sample<-rnorm(20)
 write.csv(data.frame(normal_sample),file.path(out,"normal-sample.csv"),row.names=FALSE)
 x <- solve(A,b); inverse <- rref(cbind(A,diag(10)))[,11:20]
 write.csv(A,file.path(out,"matrix.csv"),row.names=FALSE)
 write.csv(cbind(b,solution=x),file.path(out,"solution.csv"),row.names=FALSE)
 list(transpose_identity_errors=identity_errors,normal_sample_mean=mean(normal_sample),rank=qr(A)$rank,condition_number=kappa(A),solve_residual=max(abs(A%*%x-b)),inverse_residual=max(abs(A%*%inverse-diag(10))),inverse_agreement=max(abs(inverse-solve(A))),transpose_identity_error=max(abs(t(A%*%A)-t(A)%*%t(A))))
}
