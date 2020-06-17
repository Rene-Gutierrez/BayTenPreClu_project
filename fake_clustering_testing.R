d <- numeric(n)

for(i in 1:n){
  # d[i] <- sum(((y[[i]][1] - m[[1]])^2 + (y[[i]][2] - m[[1]][2])^2) - ((y[[i]][1] - m[[2]])^2 + (y[[i]][2] - m[[2]][2])^2))
  d[i] <- abs((y[[i]][1] - y[[i]][2]))
}

ord <- sort(d, index.return = TRUE)

p1 <- ord$ix[1:(n/2)]
p2 <- ord$ix[(n/2):n]


log(det(lisCov(y[p1]))) + log(det(lisCov(y[p2])))
log(det(lisCov(y[z==1]))) + log(det(lisCov(y[z==2])))

plot(matrix(unlist(y[z==1]), nrow = sum(z==1), ncol = 2, byrow = TRUE),
     col = rgb(1,0,0,0.3), pch = 16, cex = 0.05, ylim = c(-7, 7), xlim = c(-7, 7),
     xlab = "",
     ylab = "")
par(new=TRUE)
plot(matrix(unlist(y[z==2]), nrow = sum(z==2), ncol = 2, byrow = TRUE),
     col = rgb(0,0,1,0.3), pch = 16, cex = 0.05, ylim = c(-7, 7), xlim = c(-7, 7),
     xlab = "",
     ylab = "")

plot(matrix(unlist(y[p1]), nrow = length(p1), ncol = 2, byrow = TRUE),
     col = rgb(1,0,0,0.3), pch = 16, cex = 0.05, ylim = c(-7, 7), xlim = c(-7, 7),
     xlab = "",
     ylab = "")
par(new=TRUE)
plot(matrix(unlist(y[p2]), nrow = length(p2), ncol = 2, byrow = TRUE),
     col = rgb(0,0,1,0.3), pch = 16, cex = 0.05, ylim = c(-7, 7), xlim = c(-7, 7),
     xlab = "",
     ylab = "")

t(matrix(unlist(y[p2]), nrow = length(p2), ncol = 2, byrow = TRUE)) %*% matrix(unlist(y[p2]), nrow = length(p2), ncol = 2, byrow = TRUE) * 2/ n
t(matrix(unlist(y[p1]), nrow = length(p1), ncol = 2, byrow = TRUE)) %*% matrix(unlist(y[p1]), nrow = length(p1), ncol = 2, byrow = TRUE) * 2/ n

t(matrix(unlist(y[z==1]), nrow = sum(z==1), ncol = 2, byrow = TRUE)) %*% matrix(unlist(y[z==1]), nrow = sum(z==1), ncol = 2, byrow = TRUE) * 2/ n
t(matrix(unlist(y[z==2]), nrow = sum(z==2), ncol = 2, byrow = TRUE)) %*% matrix(unlist(y[z==2]), nrow = sum(z==2), ncol = 2, byrow = TRUE) * 2/ n
