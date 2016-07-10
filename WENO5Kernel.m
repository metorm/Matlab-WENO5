function [ R ] = WENO5Kernel( v1,v2,v3,v4,v5 )
%Core calculation of WENO5

	c13d12 = 13.0 / 12.0;
	c1d4 = 1.0 / 4.0;
	
	 FrontS1 = (v1 - 2.0 * v2 + v3);	 BackS1 = (v1 - 4.0 * v2 + 3.0 * v3);
	 FrontS2 = (v2 - 2.0 * v3 + v4);	 BackS2 = (v2 - v4);
	 FrontS3 = (v3 - 2.0 * v4 + v5);	 BackS3 = (3.0 * v3 - 4.0 * v4 + v5);
	 S1 = c13d12 * FrontS1 * FrontS1 + c1d4 * BackS1 * BackS1;
	 S2 = c13d12 * FrontS2 * FrontS2 + c1d4 * BackS2 * BackS2;
	 S3 = c13d12 * FrontS3 * FrontS3 + c1d4 * BackS3 * BackS3;

	 epsilon = 1e-6;

	 Alpha1 = 0.1 / ((S1 + epsilon)^2);
	 Alpha2 = 0.6 / ((S2 + epsilon)^2);
	 Alpha3 = 0.3 / ((S3 + epsilon)^2);

	 AlphaSum = Alpha1 + Alpha2 + Alpha3;

	 omega1 = Alpha1 / AlphaSum;
	 omega2 = Alpha2 / AlphaSum;
	 omega3 = Alpha3 / AlphaSum;

	 R = (omega1 / 3)*v1 + ((-7)*omega1 - omega2)*v2 / 6 ...
		+ (11 * omega1 + 5 * omega2 + 2 * omega3)*v3 / 6 ...
		+ (omega2 / 3 + 5 * omega3 / 6)*v4 + omega3*v5 / (-6);
end

