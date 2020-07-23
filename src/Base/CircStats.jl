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
    r = abs(sum(w.*exp.(im*α)))/sum(abs.(w))
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
function dirsigtest(diraxis, dirvec)
	 #  compute the trial by trial dot products (really the square of the dot product)
	dot_prods = [real(dirvec) imag(dirvec)] * [real(diraxis) imag(diraxis)]'
	 # now compute p value
	p=pvalue(OneSampleTTest(dot_prods[:], 0))

	return p
end


"""
This function calculates the ROC curve, which represents the 1-specificity and sensitivity of two classes of data, (i.e., class_1 and class_2).

The function also returns all the needed quantitative parameters: threshold position, distance to the optimum point, sensitivity,
	specificity, accuracy, area under curve (AROC), positive and negative predicted values (PPV, NPV), false negative and
	positive rates (FNR, FPR), false discovery rate (FDR), false omission rate (FOR), F1 score, Matthews correlation coefficient
	(MCC), Informedness (BM) and Markedness; as well as the number of true positives (TP), true negatives (TN), false positives (FP),
	and false negatives (FN).

Rewrite from Víctor Martínez-Cagigal (2020). ROC Curve, MATLAB Central File Exchange. Retrieved February 11, 2020.
"""
function roccurve(class_1, class_2)

	# Calculating the threshold values between the data points
	s_data = unique(sort(hcat(class_1, class_2),dims=1))  # Sorted data points
	# s_data(isnan(s_data)) = []                 # Delete NaN values
	d_data = diff(s_data)                      # Difference between consecutive points
	if isempty(d_data)
		AROC = 0
		return AROC
	    # @warn "Both class data are the same!"
	end
	d_data = vcat(d_data,d_data[end]) # Last point
	thres=[];
	thres = append!(thres,s_data[1] - d_data[1])                 # First point  # Threshold values
	thres = vcat(thres,s_data + d_data./2)
	# Calculating the sensibility and specificity of each threshold
	curve = zeros(size(thres,1),2)
	distance = zeros(size(thres,1),1)
	a=1
	for i = 1:length(thres)
	    TP = sum(class_2 .>= thres[i])    # True positives
	    FP = sum(class_1 .>= thres[i])    # False positives
	    FN = sum(class_2 .< thres[i])     # False negatives
	    TN = sum(class_1 .< thres[i])     # True negatives

	    curve[i,1] = TP/(TP + FN)   # Sensitivity
	    curve[i,2] = TN/(TN + FP)	# Specificity

	    # Distance between each point and the optimum point (0,1)
	    distance[i]= sqrt((1-curve[i,1])^2+(curve[i,2]-1)^2)
	end

	# Optimum threshold and parameters
	opt = findmin(distance)[2][1]
	TP = sum(class_2 .>= thres[opt])    # No. true positives
	FP = sum(class_1 .>= thres[opt])    # No. false positives
	FN = sum(class_2 .< thres[opt])     # No. false negatives
	TN = sum(class_1 .< thres[opt])     # No. true negatives

	# Output parameters
	# param = DataFrame()
	# param.Threshold = thres[opt]                #Optimum threshold position
	# param.Sensi = curve[opt,1]                 # Sensitivity
	# param.Speci = curve[opt,2]                 # Specificity
	# param.AROC  = abs(trapz(1 .- curve[:,2], curve[:,1])) # Area under curve
	# param.Accuracy = (TP+TN)/(TP+TN+FP+FN)     # Aaccuracy
	# param.PPV   = TP/(TP+FP)                   # Positive predictive value
	# param.NPV   = TN/(TN+FN)                   # Negative predictive value
	# param.FNR   = FN/(FN+TP)                   # False negative rate
	# param.FPR   = FP/(FP+TN)                   # False positive rate
	# param.FDR   = FP/(FP+TP)                   # False discovery rate
	# param.FOR   = FN/(FN+TN)                   # False omission rate
	# param.F1_score = 2*TP/(2*TP+FP+FN)         # F1 score
	# param.MCC   = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))  # Matthews correlation coefficient
	# param.BM    = curve[opt,1] + curve[opt,2] - 1    # Informedness
	# param.MK    = TP/(TP+FP) + TN/(TN+FN) - 1        # Markedness
	#
	# param.TP = TP    # No. true positives
	# param.FP = FP    # No. false positives
	# param.FN = FN    # No. false negatives
	# param.TN = TN    # No. true negatives
    AROC  = abs(trapz(1 .- curve[:,2], curve[:,1]))  # so far , only need AUC
	return AROC
end
