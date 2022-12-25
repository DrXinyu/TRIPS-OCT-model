function Stokes = tom2Stokes(Hc,Vc)
   
    S1 = abs(Hc).^2 - abs(Vc).^2;
    S2 = 2*real(Hc.*conj(Vc));
    S3 = -2*imag(Hc.*conj(Vc));
    dimtom = size(Hc);
    Stokes = squeeze(reshape(cat(4,S1,S2,S3),[dimtom,3]));
    
end

