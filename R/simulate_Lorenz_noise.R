#' Simulate the Lorenz System with noise
#' @import stats
#' @import deSolve
#' @param N Length of Times Series
#' @param delta_t Step size for time scale. If \code{NULL} this is derived using \code{N} and \code{tmax}. \code{DEFAULT} to \code{NULL}.
#' @param tmax Upper bound of the time scale. This argument is ignored if \code{delta_t} is provided. \code{DEFAULT} to \code{50}.
#' @param X0 Initial value for X at t=0. \code{DEFAULT} to 0.
#' @param Y0 Initial value for Y at t=0. \code{DEFAULT} to 1.
#' @param Z0 Initial value for Z at t=0. \code{DEFAULT} to 1.
#' @param sdX Use this argument to rescale the X-coordinate to have the desired standard deviation (exactly). This is ignored if set to \code{NULL}. \code{DEFAULT} to \code{NULL}.
#' @param sdY Use this argument to rescale the Y-coordinate to have the desired standard deviation (exactly). This is ignored if set to \code{NULL}. \code{DEFAULT} to \code{NULL}.
#' @param sdZ Use this argument to rescale the Z-coordinate to have the desired standard deviation (exactly). This is ignored if set to \code{NULL}. \code{DEFAULT} to \code{NULL}.
#' @param sdnoiseX Standard deviation of Gaussian noise of X-coordinate. If set to \code{0}, no noise is created.
#' @param sdnoiseY Standard deviation of Gaussian noise of Y-coordinate. If set to \code{0}, no noise is created.
#' @param sdnoiseZ Standard deviation of Gaussian noise of Z-coordinate. If set to \code{0}, no noise is created.
#' @param s s-parameter of the Lorenz ODE. See Vignette for further details. \code{DEFAULT} to \code{10}, which is the original value chosen by Lorenz.
#' @param r r-parameter of the Lorenz ODE. See Vignette for further details. \code{DEFAULT} to \code{28}, which is the original value chosen by Lorenz.
#' @param b b-parameter of the Lorenz ODE. See Vignette for further details. \code{DEFAULT} to \code{8/3}, which is the original value chosen by Lorenz.
#' @param naive Logical whether naive calculation should be used. \code{DEFAULT} to \code{FALSE}.
#' @param return_time Logical whether the time-coordinate should be included in the returned \code{data.frame}. \code{DEFAULT} to \code{TRUE}.
#' @returns Returns a three dimensional time series as \code{data.frame} following the Lorenz system (Lorenz, 1963, <doi:10.1175/1520-0469(1963)020<0130:DNF>2.0.CO;2>).
#' @references Lorenz, E. N. (1963). Deterministic nonperiodic flow. Journal of atmospheric sciences, 20(2), 130-141. <doi:10.1175/1520-0469(1963)020<0130:DNF>2.0.CO;2>
#' @export

simulate_Lorenz_noise <- function(N = 1000, delta_t = NULL, tmax = 50,
                         X0 = 0, Y0 = 1, Z0 = 1,
                         sdX = NULL, sdY = NULL, sdZ = NULL,
                         sdnoiseX, sdnoiseY, sdnoiseZ,
                         s = 10, r = 28, b = 8/3,
                         naive = FALSE,
                         return_time = TRUE)
{

     ### prepare function ----------------------------------------------------
     if(!is.null(delta_t)){
          if(delta_t <= 0) stop("delta_t needs to be numeric and positive.")
          tmax <- (N-1)*delta_t
     }
     times <- seq(0, tmax, length.out = N)

     state <- c("X" = X0, "Y" = Y0, "Z" = Z0)
     parameters <- c("s" = s, "r" = r, "b" = b)

     Lorenz <- function(t, state, parameters) {
          with(as.list(c(state, parameters)), {
               dX <- s * (Y - X)
               dY <- X * (r - Z) - Y
               dZ <- X * Y - b * Z
               list(c(dX, dY, dZ))
          })
     }

     ### solve ODE ----------------------------------------------------------------
     if(naive)
     {

          NaiveLorenz <- function(times, state, parameters)
          {
               LorenzODE <- function(dt, X, Y, Z, s, r, b){
                    dX <- s * (Y - X)
                    dY <- X * (r - Z) - Y
                    dZ <- X * Y - b * Z
                    return(c(dX*dt + X, dY*dt + Y, dZ*dt + Z))
               }

               XYZ <- rbind(c(X0, Y0, Z0), matrix(NA, nrow = length(times)-1, ncol = 3))
               for(i in 2:length(times))
               {
                    dt <- times[i]-times[i-1]
                    XYZ[i,] <- LorenzODE(dt = dt, X = XYZ[i-1,1], Y = XYZ[i-1,2], Z = XYZ[i-1,3],
                                      s = s, r = r, b = b)
               }
               XYZ <- data.frame(times, XYZ)
               names(XYZ) <- c("time", "X", "Y", "Z")
               return(XYZ)
          }
          out <- NaiveLorenz(times = times, state = state, parameters = parameters)
     }else{
          out <- deSolve::ode(y = state, times = times, func = Lorenz, parms = parameters)  |> data.frame()
     }

     if(all(out[,-1] == 0)) stop("The point (0,0,0) as a starting point results in a constant (=0) time series with noise.")

     ### rescale results ----------------------------------------------------------
     if(!is.null(sdX))
     {
          if(!is.numeric(sdX) | sdX < 0) stop("sdX needs to be positive and numeric.")
          out$X <- scale(out$X) * sdX
     }
     if(!is.null(sdY))
     {
          if(!is.numeric(sdY) | sdY < 0) stop("sdY needs to be positive and numeric.")
          out$Y <- scale(out$Y) * sdY
     }
     if(!is.null(sdZ))
     {
          if(!is.numeric(sdZ) | sdZ < 0) stop("sdZ needs to be positive and numeric.")
          out$Z <- scale(out$Z) * sdZ
     }

     ### add noise ---------------------------------------------------------------
     Noise <- mvtnorm::rmvnorm(n = N, mean = c(0, 0, 0),
                               sigma = diag(c(sdnoiseX^2, sdnoiseY^2, sdnoiseZ^2)))
     out[,-1] <- out[,-1] + Noise

     if(!return_time) out$time <- NULL
     return(out)
}


