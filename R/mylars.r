mylars<-function (x, y, l1.weights=NULL, design=TRUE, 
    type = c("lasso", "lar", "forward.stagewise"), 
    trace = FALSE, Gram=NULL, eps = .Machine$double.eps,
    max.steps=NULL, use.Gram = TRUE) 
{
    call <- match.call()
    type <- match.arg(type)
    TYPE <- switch(type, lasso = "LASSO", lar = "LAR", forward.stagewise = "Forward Stagewise")
    if (trace) 
        cat(paste(TYPE, "sequence\n"))
    nm <- dim(x)
    n <- nm[1]
    m <- nm[2]
    im <- inactive <- seq(m)
    one <- rep(1, n)
    vn <- dimnames(x)[[2]]
    meanx <- drop(one %*% x)/n
    if (is.null(l1.weights)==FALSE) x<-scale(x,FALSE,l1.weights)
#    x <- scale(x, meanx, FALSE)
    normx <- sqrt(drop(one %*% (x^2)))
    nosignal <- normx/sqrt(n) < eps
    if (any(nosignal)) {
        ignores <- im[nosignal]
        inactive <- im[-ignores]
        normx[nosignal] <- eps * sqrt(n)
        if (trace) 
            cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < eps; dropped for good\n")
    }
    else ignores <- NULL
    names(normx) <- NULL
#   x <- scale(x, FALSE, normx)
    if (use.Gram & (is.null(Gram)==TRUE)) {
        if (m > 500 && n < m) 
            cat("There are more than 500 variables and n<m;\nYou may wish to restart and set use.Gram=FALSE\n")
        if (trace) 
            cat("Computing X'X .....\n")
        Gram <- t(x) %*% x
	#print(Gram)
    }
    mu <- mean(y)
    #y <- drop(y - mu)
    y<-drop(y); 
    Cvec <- drop(t(y) %*% x)
    #print(Cvec); 
    ssy <- sum(y^2)
    residuals <- y
    if (is.null(max.steps)) 
        max.steps <- 8 * min(m, n - 1)
    beta <- matrix(0, max.steps + 1, m)
    Gamrat <- NULL
    arc.length <- NULL
    R2 <- 1
    RSS <- ssy
    first.in <- integer(m)
    active <- NULL
    actions <- as.list(seq(max.steps))
    drops <- FALSE
    Sign <- NULL
    R <- NULL
    k <- 0
    while ((k < max.steps) & (length(active) < min(m - length(ignores), 
        n - 1))) {
        action <- NULL
        k <- k + 1
        C <- Cvec[inactive]
        Cmax <- max(abs(C))
        if (!any(drops)) {
            new <- abs(C) >= Cmax - eps
            C <- C[!new]
            new <- inactive[new]
            for (inew in new) {
                if (use.Gram) {
                  R <- updateR(Gram[inew, inew], R, drop(Gram[inew, 
                    active]), Gram = TRUE, eps = eps)
                }
                else {
                  R <- updateR(x[, inew], R, x[, active], Gram = FALSE, 
                    eps = eps)
                }
                if (attr(R, "rank") == length(active)) {
                  nR <- seq(length(active))
                  R <- R[nR, nR, drop = FALSE]
                  attr(R, "rank") <- length(active)
                  ignores <- c(ignores, inew)
                  action <- c(action, -inew)
                  if (trace) 
                    cat("LARS Step", k, ":\t Variable", inew, 
                      "\tcollinear; dropped for good\n")
                }
                else {
                  if (first.in[inew] == 0) 
                    first.in[inew] <- k
                  active <- c(active, inew)
                  Sign <- c(Sign, sign(Cvec[inew]))
                  action <- c(action, inew)
                  if (trace) 
                    cat("LARS Step", k, ":\t Variable", inew, 
                      "\tadded\n")
                }
            }
        }
        else action <- -dropid
        Gi1 <- backsolve(R, backsolvet(R, Sign))
        dropouts <- NULL
        if (type == "forward.stagewise") {
            directions <- Gi1 * Sign
            if (!all(directions > 0)) {
                if (use.Gram) {
                  nnls.object <- nnls.lars(active, Sign, R, directions, 
                    Gram[active, active], trace = trace, use.Gram = TRUE, 
                    eps = eps)
                }
                else {
                  nnls.object <- nnls.lars(active, Sign, R, directions, 
                    x[, active], trace = trace, use.Gram = FALSE, 
                    eps = eps)
                }
                positive <- nnls.object$positive
                dropouts <- active[-positive]
                action <- c(action, -dropouts)
                active <- nnls.object$active
                Sign <- Sign[positive]
                Gi1 <- nnls.object$beta[positive] * Sign
                R <- nnls.object$R
                C <- Cvec[-c(active, ignores)]
            }
        }
        A <- 1/sqrt(sum(Gi1 * Sign))
        w <- A * Gi1
        if (!use.Gram) 
            u <- drop(x[, active, drop = FALSE] %*% w)
        if (length(active) >= min(n - 1, m - length(ignores))) {
            gamhat <- Cmax/A
        }
        else {
            if (use.Gram) {
                a <- drop(w %*% Gram[active, -c(active, ignores), 
                  drop = FALSE])
            }
            else {
                a <- drop(u %*% x[, -c(active, ignores), drop = FALSE])
            }
            gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
            gamhat <- min(gam[gam > eps], Cmax/A)
        }
        if (type == "lasso") {
            dropid <- NULL
            b1 <- beta[k, active]
            z1 <- -b1/w
            zmin <- min(z1[z1 > eps], gamhat)
            if (zmin < gamhat) {
                gamhat <- zmin
                drops <- z1 == zmin
            }
            else drops <- FALSE
        }
        beta[k + 1, ] <- beta[k, ]
        beta[k + 1, active] <- beta[k + 1, active] + gamhat * 
            w
        if (use.Gram) {
            Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*% 
                w
        }
        else {
            residuals <- residuals - gamhat * u
            Cvec <- drop(t(residuals) %*% x)
        }
        Gamrat <- c(Gamrat, gamhat/(Cmax/A))
        arc.length <- c(arc.length, gamhat)
        if (type == "lasso" && any(drops)) {
            dropid <- seq(drops)[drops]
            for (id in rev(dropid)) {
                if (trace) 
                  cat("Lasso Step", k + 1, ":\t Variable", active[id], 
                    "\tdropped\n")
                R <- downdateR(R, id)
            }
            dropid <- active[drops]
            beta[k + 1, dropid] <- 0
            active <- active[!drops]
            Sign <- Sign[!drops]
        }
        if (!is.null(vn)) 
            names(action) <- vn[abs(action)]
        actions[[k]] <- action
        inactive <- im[-c(active, ignores)]
    }
    beta <- beta[seq(k + 1), ]
    dimnames(beta) <- list(paste(0:k), vn)
    if (is.null(l1.weights)==FALSE) beta<-scale(beta,FALSE,l1.weights)
if (design==TRUE) {
    if (trace) 
        cat("Computing residuals, RSS etc .....\n")
    residuals <- y - x %*% t(beta)
 #  beta <- scale(beta, FALSE, normx)
    RSS <- apply(residuals^2, 2, sum)
    R2 <- 1 - RSS/RSS[1]
    Cp <- ((n - k - 1) * RSS)/rev(RSS)[1] - n + 2 * seq(k + 1)
    } else {RSS<-R2<-Cp<-NULL}
    object <- list(call = call, type = TYPE, R2 = R2, RSS = RSS, 
        Cp = Cp, actions = actions[seq(k)], entry = first.in, 
        Gamrat = Gamrat, arc.length = arc.length, Gram = if (use.Gram) Gram else NULL, 
        beta = beta, mu = mu, normx = normx, meanx =
	meanx,l1.weights=l1.weights)
    class(object) <- "lars"
    object
}
