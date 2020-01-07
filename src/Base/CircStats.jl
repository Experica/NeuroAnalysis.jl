"Mean for circular data"
circmean(α::AbstractVector,w=ones(length(α)))=sum(w.*exp.(im*α))
"""
Mean Resultant Vector Length for circular data
    α: sample of angles in radians
    w: number of incidences in case of binned angle data
    d: spacing of bin centers for binned data, if supplied
       correction factor is used to correct for bias
"""
function circr(α::AbstractVector,w=ones(length(α)),d=nothing)
    r = abs(sum(w.*exp.(im*α)))/sum(w)
    if d!=nothing
        r *= d/2/sin(d/2)
    end
    r
end

"Circular Variance"
circvar(α::AbstractVector,w=ones(length(α)),d=nothing) = 1-circr(α,w,d)

"""
Transforms p-axial data to a common scale
  α: sample of angles in radians
  p: number of modes
"""
circaxial(α,p=2)=(p*α)%(2pi)

"Mean direction for circular data with axial"
circaxialmean(α::AbstractVector,p=2)=angle(sum(exp.(im*α*p)))/p

"""
Hotelling's T-Squared test for one multivariate sample.

Performs Hotelling's T^2 test on multivariate samples X to determine
if the data have mean MU.  X should be a NxP matrix with N observations
of P-dimensional data, and the mean MU to be tested should be 1xP.
the significance level α is 0.05 by default.

H is 1 if the null hypothesis (that the mean of X is equal to MU) can be
rejected at significance level α.  P is the actual P value.

The code is based on HotellingT2.m by A. Trujillo-Ortiz and
R. Hernandez-Walls.
"""
function hotellingt2test(X,mu,α = 0.05)
	n=size(X,1); p=size(X,2); m = mean(X,dims=1); S = cov(X);
	T2 = n*(m-mu)*inv(S)*(m-mu)'  # T-square statistics

	if n>= 50 # Chi-square approximation
		P = 1-chisqcdf(T2[1],p)
	else  # F approximation
		F = (n-p)/((n-1)*p)*T2[1]
		P = ccdf(FDist(p, n-p), F)
	end
	H = P<α
	return P
end

"""
Direction tuning significance using dot product with empirical orientation preference
This function calculates the probability that the "true" direction tuning
vector of a neuron has non-zero length. It performs this by empirically
determing the unit orientation vector, and then computing the dot
product of the direction vector for each trial onto the overall orientation
vector, and then looking to see if the average is non-zero.

Inputs:  ANGLES is a vector of direction angles at which the response has
       been measured.
       RATES is the response of the neuron in to each angle; each row
       should contain responses from a different trial.
Output:  P the probability that the "true" direction tuning vector is
       non-zero.

See: Mazurek, Kagan, Van Hooser 2014;  Frontiers in Neural Circuits
"""
function dirsigtest(orivec, dirvec)
	 #  compute the trial by trial dot products (really the square of the dot product)
	dot_prods = [real(dirvec) imag(dirvec)] * [real(orivec) imag(orivec)]'
	 # now compute p value
	p=pvalue(OneSampleTTest(dot_prods[:], 0))

	return p
end
