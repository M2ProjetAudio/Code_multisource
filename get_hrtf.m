function [hrir,H,P,V]=get_hrtf(Lframe,B,az,Nbfreq,Ntheta,freqIndexes)


hrtffile = xml.dbGetFile(...
    'impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa');
hrir_SOFA= SOFAload(hrtffile);

HRIR = hrir_SOFA.Data.IR;
HRIR = permute(HRIR,[3 2 1]);

HRTF = NaN(Lframe,2,360);
% compute the anechoic hrtfs by fourier transforming irs. zero padd irs if
% irs length is smaller than a frame length, otherwise truncate irs.

if Lframe <= 2048
    for i=1:length(az)
        HRTF(:,:,i)=fft(HRIR(1:Lframe,:,i));
    end
else
    for i=1:length(az)
        HRTF(:,:,i)=fft([HRIR(:,:,i) ; zeros(Lframe-2048,2)]);
    end
end

H = squeeze(HRTF(:,2,:)./HRTF(:,1,:));

P=NaN(2,2,Nbfreq,Ntheta);

for k=1:B
    for ntheta=1:Ntheta
        V=[1;H(freqIndexes(k),ntheta)];
        P(:,:,k,ntheta)=(V*V')/(V'*V); %donnee
    end
end
hrir = simulator.DirectionalIR(xml.dbGetFile...
    ('impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa'));


