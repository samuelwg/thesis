function [minimumPhaseSignal,realCepstrum] = homoDeco(s) 
% [sm] = mps(s) 
% create minimum-phase spectrum sm from complex spectrum s 
%
realCepstrum= real(ifft(log(abs(fft(s)))));
N=length(realCepstrum);

[sx,sy]=size(realCepstrum);

filtre=2*ones(sx,sy);
filtre(1)=1;

if mod(N,2)==0
	filtre(  N/2+1 )=1;	
	filtre( (N/2+2):N)=0;	
else
	filtre( ((N+1)/2+1):N ) =0;
endif

complexCeptrum = filtre.*realCepstrum;

minimumPhaseSignal = real(ifft(exp(fft(complexCeptrum))));

