# ==========================================================================
# Author: Luke Lloyd-Jones
# Date started: 10/10/2016
# Date last updated: 09/02/2018
# Function takes summary statistics from a linear model or
# linear mixed model genome-wide association study and
# maps the genetic effect estimates to the odds ratio. The
# mapped results are appended to the original data frame.
# ==========================================================================
LmToOddsRatio <- function(lmm.fil, k, std.err)
{
  # ------------------------------------------------------------------------
  # Args:
  #   lmm.fil: Data frame with genetic effect (BETA), 
  #            allele frequency (FREQ) and sample prevalence (K)
  #         k: If the user just want to specify one k for all 
  #            genetic variants
  #   std.err: Default should be 0. Logical for if the standard error 
  #            transformation is required. std.err = 1 implies the SE 
  #            transformation will be returned. WARNING: REQUIRES A 
  #            "SE" COLUMN AND N = NUMBER OF INDIVIDUALS COLUMN
  # 
  # Returns:
  #   Original lmm.fil data frame with the odds ratio appended
  # ------------------------------------------------------------------------
  # Grab the columns names
  df.names <- colnames(lmm.fil) 
  # Get the column name for the allele frequency 
  p.col    <- which(df.names == "Freq" | df.names == "FREQ" | 
                    df.names == "FRQ"  |  df.names == "A1FREQ" | 
                    df.names == "af"| df.names == "freq")
  # If one of the column names is not in the set of possible names request 
  # the user to change the column name
  if (length(p.col) == 0)
  {
    stop("Please rename the allele frequency column FREQ or add an allele 
         frequency for the effect allele from a reference. Alternatively, 
         check that you have the correct seperator for your results file.")
  }
  # --------------------------------------------------- 
  # Get the column name for the reference allele effect
  # --------------------------------------------------- 
  p.beta <- which(df.names == "b" | df.names == "BETA" | df.names == "beta" | 
                  df.names == "EFFECT")
  # If one of the column names is not in the set of possible names 
  # request the user to change the column name
  if (length(p.beta) == 0)
  {
    stop("Please rename the genetic effect column BETA or check that you 
         have the correct seperator for your results file")
  }
  if (length(p.beta) > 1)
  {
    stop("You have multiple genetic coefficent columns please just submit one")
  }
  # Check to see if the effects have any NAs. Remove and notify user
  if (length(which(is.na(lmm.fil[, p.beta]))) != 0)
  {
    print("Removing results with NAs for the genetic effect")
    lm.fil <- lmm.fil[-which(is.na(lmm.fil[, p.beta])), ]
  }
  # Check to see if the allele frequencies have any NAs. 
  # Remove and notify user
  if (length(which(is.na(lmm.fil[, p.col]))) != 0)
  {
    print("Removing results with no allele frequency")
    lm.fil <- lmm.fil[-which(is.na(lmm.fil[, p.beta])), ]
  }
  # Set the variables for the allele and the effect
  p    <- lmm.fil[, p.col]
  beta <- lmm.fil[, p.beta]
  # ------------------------
  # Odds ratio 3. A remnant.
  # ------------------------
  p_1 <- p  + (beta * p * (1 - p)) / k
  p_0 <- (p - k * p_1) / (1 - k)
  or.map.1  <- (p_1 * (1 - p_0)) / (p_0 * (1 - p_1))
  # -----------------------
  # Solve for the quadratic
  # -----------------------
  pp <- p
  c1 <- (1 - k) * (3 * beta - 2 * k * beta) + (beta * ((1 - k) ^ 2) * 
        (1 - 2 * k)) / k
  c2 <- (1 - k) * (1 - 4 * beta * pp) - (2 * beta * pp * (1 - k) * 
        (1 - 2 * k)) / k
  c3 <- (pp ^ 2) * beta * (1 - 2 * k) / k + pp * (k - 1 + beta)
  f0 <- c1 + c2 + c3
  # ------------------------------------------------------------------------
  # Calculate the upper and lower bounds for returning sensible ORs
  # ------------------------------------------------------------------------
  lb <- (p * k - p * (k ^ 2) - k * (1 - k)) / 
        ((1 - k) - 2 * p * (1 - k) + p ^ 2 - 2 * (p ^ 2) * k + p * k)
  ub <- (p * k - p * (k ^ 2)) / (p ^ 2 - 2 * (p ^ 2) * k + p * k)
  # Set up an empty array
  or.map.gld <- array(0, length(p))
  if (length(which(is.na(c3))) != 0 & length(which(is.na(f0))) != 0)
  {
    # --------------------
    # Odds ratio quadratic
    # --------------------  
    p.0.1 <- (-c2 + sqrt(c2 ^ 2 - 4 * c1 * c3)) / (2 * c1)
    p.0.2 <- (-c2 - sqrt(c2 ^ 2 - 4 * c1 * c3)) / (2 * c1)
    # There should be only one solution in (0, 1). Take the set that 
    # has all elements in (0, 1)
    if (all(p.0.1 > 0 & p.0.1 < 1))
    {
      p_0 <- p.0.1
    } else if (all(p.0.2 > 0 & p.0.2 < 1))
    {
      p_0 <- p.0.2
    }
    p_1 <- (p - (1 - k) *  p_0 ) / k  
    or.map.gld  <- (p_1 * (1 - p_0)) / (p_0 * (1 - p_1))
    or.map.gld[which(is.na(c3))] <- NA
    or.map.gld[which(is.na(f0))]  <- NA   
  } else if (!all(beta > lb & beta < ub))  
  {
    c1[which(beta < lb | beta> ub)] <- NA
    c2[which(beta < lb | beta> ub)] <- NA
    c3[which(beta < lb | beta> ub)] <- NA
    p.0.1 <- (-c2 + sqrt(c2 ^ 2 - 4 * c1 * c3)) / (2 * c1)
    p.0.2 <- (-c2 - sqrt(c2 ^ 2 - 4 * c1 * c3)) / (2 * c1)
    # There should be only one solution in (0, 1). 
    # Take the set that has all elements in (0, 1)
    if (!all(is.na(p.0.1)) & !all(is.na(p.0.1)))
    {
      if (all(p.0.1[-which(is.na(p.0.1))] > 0 & 
              p.0.1[-which(is.na(p.0.1))] < 1))
      {
        p_0 <- p.0.1
      } else if (all(p.0.2[-which(is.na(p.0.2))] > 0 & 
                     p.0.2[-which(is.na(p.0.2))] < 1))
      {
        p_0 <- p.0.2
      }
    } else
    {
      p_0 <- NA
    }
    p_1 <- (p - (1 - k) *  p_0 ) / k
    or.map.gld  <- (p_1 * (1 - p_0)) / (p_0 * (1 - p_1))
  } else
  {
    # If there are no concerns with NAs and ORs on the boundary
    # calculate the full set of results
    or.map.gld  <- (p_1 * (1 - p_0)) / (p_0 * (1 - p_1))
  }
  # ------------------------------------------------------------------------
  # If standard error mappings are requested
  # ------------------------------------------------------------------------
  if (std.err == 1)
  {
    p.se <- which(df.names == "SE" | df.names == "se")
    p.n  <- which(df.names == "N" | df.names == "NMISS" | 
                  df.names == "nmiss")
    # If one of the column name is not in the set of possible names request 
    # the user to change the column name
    if (length(p.se) == 0)
    {
      stop("Please rename the standard error of the genetic effect column 
           SE or check that you have the correct seperator for your results 
           file")
    }
    if (length(p.n) == 0)
    {
      stop("Please rename the number of individuals column N or check 
           that you have the correct seperator for your results file")
    }
    se <- lmm.fil[, p.se]
    n  <- lmm.fil[, p.n]
    d  <- (beta ) / (2 * ((se ^ 2) * (n - 2) + (beta ^ 2)))
    a  <- -beta
    b  <- beta - 2 * d * k * beta
    c  <- beta * 2 * k * (1 - k) * (d ^ 2) - k * (1 - k) * d + 
          k * beta * (d - (d^2))
    p0_est.1 <- (-b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a)
    p0_est.2 <- (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a)  
    p0.prop  <- array(NA, length(p0_est.1))
    p0.prop[which(beta < 0)]   <- p0_est.1[which(beta < 0)]
    p0.prop[which(beta >= 0)]  <- p0_est.2[which(beta >= 0)]
    p1_est <- p0.prop + d
    or.est.2 <- (p1_est - p1_est * p0.prop) / (p0.prop - p1_est * p0.prop)
    # Calculate the SE quadratic to evaluate when the transformation 
    # will break down
    var <- se^2
    a   <- ((n - 2) ^ 2)
    b   <- (2 * (beta ^ 2) * (n - 2) - 2 * k * (n - 2) + 2 * (k^2) * (n-2))
    c   <- (beta^4) + (beta^2) * (k^2) - (beta^2) * k
    sd  <- var^2 * a + var * b + c
    if (!all(sd > 0))
    {
      or.est.2[which(sd < 0)] <- NA
    }
    # Return the results with the SE mapping
    res <- data.frame(lmm.fil,  or.est.2)
    colnames(res) <- c(colnames(lmm.fil), "OR_SE")
    return(res)
    }
  # -------------------
  # Bind up the results 
  # -------------------
  if (std.err != 1)
  {
    res <- data.frame(lmm.fil, or.map.gld)
    colnames(res) <- c(colnames(lmm.fil), "OR")
    return(res)
  }
}




