more off;
order = 6; 

% First run the first if to generate the filters only once, then the second
if 0
	t0=-1;
	tfinal = 40;
	dt = 0.0001;
	nPoints = round((tfinal-t0)/dt)
	tRange = t0+(1:nPoints)*dt;

	tRange1 = tRange( find(tRange<0)  );
	tRange2temp = tRange( find(tRange >= 0));
	tRange2 = tRange2temp( find(tRange2temp<=pi/2)  );
	htotal = zeros(order+1,nPoints);
	for l=0:6
		l
		h1 = legendre(l,tRange1)(1,:);
		h2 = cos(tRange2).*legendre(l,sin(tRange2))(1,:);
		h = [h1 h2 zeros(1,nPoints-length(tRange1)-length(tRange2))];
		htotal(l+1,:)=h;
	endfor
	save "h.txt" htotal;
endif

if 1
	for l=0:order
		load "h.txt";
		h = htotal(l+1,:);
		clear htotal;
		maxim = max(abs( h - homoDeco(h) ))
		if 0
			figure;
			hold on;
			plot(h(1:round(length(h)/2)),"b")		
			plot(homoDeco(h)(1:round(length(h)/2)),"r")
			hold off;
		endif		
	endfor
endif
