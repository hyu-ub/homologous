#' Test of homologous hypothesis
#'
#' Perform the test of homologous hypothesis for comparing growth curves
#'
#' @details Perform the test of homologous hypothesis for comparing growth curves
#'
#' @param x a data.frame with four columns: x, the measurement such as tumor volume; time, the time of measurement; subject, the subject id; group, the subject's treatment group
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#'
#' @return object with the local models incorporated
#'
#' @author Han Yu
#'
#' @export

homologous_test <- function(x, alternative = c("two.sided", "less", "greater")){

  if (!alternative %in% c("two.sided", "less", "greater")) {
    warning("alternative is not one of two.sided, less, or greater, using the default option of two.sided")
    alternative <- "two.sided"
  }
  
  df <- x
  time_list <- sort(unique(df$time))
  groups <- sort(unique(df$group))
  df$id <- apply(df[, c("group", "subject")], 1, paste0, collapse = "-")

  s <- d <- t_value <- df_value <- p_value <- c()

  for(i in 2:length(time_list)) {

    y0 <- df$x[df$time == time_list[i-1]]
    g0 <- df$group[df$time == time_list[i-1]]
    s0 <- df$id[df$time == time_list[i-1]]

    df_0 <- data.frame(id = s0, y0, group = g0)

    y1 <- df$x[df$time == time_list[i]]
    g1 <- df$group[df$time == time_list[i]]
    s1 <- df$id[df$time == time_list[i]]

    df_1 <- data.frame(id = s1, y1)

    df_a <- merge(df_0, df_1, by = "id")
    df_a$group <- as.numeric(as.factor(df_a$group))-1
    df_a$intercept <- 1

    df_a$itr <- df_a$y0 * df_a$group

    X <- as.matrix(df_a[,c("intercept", "group", "y0", "itr")])
    Y <- df_a[, "y1"]
    n <- nrow(X)

    y00 <- df_a$y0[df_a$group==0]
    y01 <- df_a$y0[df_a$group==1]
    y10 <- Y[df_a$group==0]
    y11 <- Y[df_a$group==1]

    if(var(y00, na.rm=TRUE)==0 | var(y01, na.rm=TRUE)==0 | var(y10, na.rm=TRUE)==0 | var(y11, na.rm=TRUE)==0) {
      s[i] <- NA
      d[i] <- NA
      t_value[i] <- NA
      df_value[i] <- NA
      p_value[i] <- NA
      next
    }

    beta <- solve(t(X) %*% X) %*% t(X) %*% Y
    # lim_fit <- lm(y1 ~ y0+group+itr, df_a)
    # summary(lim_fit)

    z <- c(0, 1, mean(y01) - mean(y00), mean(y01))

    MSE <- (sum(Y^2) - t(beta)%*%t(X)%*%Y)/(n-4)
    s2_beta <- MSE[1,1]*solve(t(X)%*%X)
    s2 <- t(z) %*% s2_beta %*% z

    s[i] <- sqrt(s2)[1,1]
    d[i] <- sum(z*beta)
    t_value[i] <- d[i]/s[i]
    df_value[i] <- n-4



    if(alternative == "two.sided") p_value[i] <- 2*pt(-abs(t_value[i]), df_value[i])
    if(alternative == "less") p_value[i] <- pt(t_value[i], df_value[i])
    if(alternative == "greater") p_value[i] <- pt(-t_value[i], df_value[i])
  }

  df_out <- data.frame(time = time_list, Estimate = d, SE = s, t_value, DF = df_value, p_value)
  return(df_out)
}
