% computes the pdf for the univariate location scale t distribution
function yvals = locsctpdf( x, df, mu, sigma )

normconst = @( df, sigma ) gamma( (df+1)/2 )/( sigma*sqrt(df*pi)*gamma(df/2) );

kernel = @( x, df, mu, sigma ) ( 1 + ( x - mu ).^2 ./ (df*sigma^2) ).^(-(df+1)/2);

yvals = normconst(df,sigma)*kernel(x,df,mu,sigma);
