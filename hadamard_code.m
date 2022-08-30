clearvars; close all; clc;

nBits = 999; % number of bits
nOfSim = 10; % number of monte carlo simulations
snr = -5:2.5:10; %signal to noise ratios
%8 bits assigned for each 3 bits
code_rate = 3/8; %code rate
H_4 = [1 1 1 1;1 -1 1 -1;1 1 -1 -1;1 -1 -1 1];
H_8 = [H_4 H_4;H_4 -H_4];%hadamard matrix of order 8
three_bit_matrix = [0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1];

%%initializations
decoded_hadamard_index = zeros(1,nBits/3);
index = zeros(1,nBits/3);
soft_decoded_hadamard_data = zeros(1,nBits);
BER = zeros(1,nOfSim);
BER_hadamard_hard = zeros(1,nOfSim);
BER_hadamard_soft = zeros(1,nOfSim);
BER_SNR = zeros(1,length(snr));
BER_SNR_hadamard_hard = zeros(1,length(snr));
BER_SNR_hadamard_soft = zeros(1,length(snr));


decoded_signal = zeros(1,nBits);
decoded_signal_hadamard = zeros(1,nBits);

for i = 1:length(snr)
    %%monte carlo simulations
    for j = 1:nOfSim

        data = randi([0 1],1,nBits); % bit stream to send
        tx = real(pskmod(data,2,pi)); %binary psk modulated signal
        data_3tuple = reshape(data,[3,nBits/3])';
        dec_data_3tuple = bi2de(data_3tuple,'left-msb')';
        tx_hadamard = H_8(dec_data_3tuple(1:end)+1,:);
        tx_hadamard= reshape(tx_hadamard',[1, nBits*8/3]);
        rx = awgn(tx,snr(i),'measured');
        rx_hadamard = awgn(tx_hadamard,snr(i)+log10(code_rate),'measured');
       %%hard decoding 
        for k= 1:nBits                  
            if ((rx(k)<0))
                decoded_signal(k) = 0;
            else
                decoded_signal(k) = 1;
            end
        end
        for k = 1:nBits*8/3
            if ((rx_hadamard(k)<0))
                decoded_signal_hadamard(k) = -1;
            else
                decoded_signal_hadamard(k) = 1;
            end
        end

        decoded_signal_hadamard = reshape(decoded_signal_hadamard,[8, nBits/3])';
        for k = 1:nBits/3
            [~,decoded_hadamard_index(k)] = min(sum(abs(decoded_signal_hadamard(k,:) - H_8)/2,2));
        end
        harddecoded_hadamard_data = three_bit_matrix(decoded_hadamard_index(1:end),:);
        harddecoded_hadamard_data = reshape(harddecoded_hadamard_data',[1,nBits]);
        %%soft decoding
        rx_hadamard = reshape(rx_hadamard,[8, nBits/3])';
        for k = 1:nBits/3
            [~,index(k)] = max(sum(rx_hadamard(k,:).*H_8,2));
        end
        soft_decoded_hadamard_data = three_bit_matrix(index(1:end),:);
        soft_decoded_hadamard_data = reshape(soft_decoded_hadamard_data',[1 nBits]);
        %%BER for each simulation
        BER(j) = sum(mod(data + decoded_signal,2))/nBits;
        BER_hadamard_hard(j) = mean(mod(data + harddecoded_hadamard_data,2));
        BER_hadamard_soft(j) = mean(mod(data + soft_decoded_hadamard_data,2));
    end
    %%average BERs
    BER_SNR(i) = mean(BER);
    BER_SNR_hadamard_hard(i) = mean(BER_hadamard_hard);
    BER_SNR_hadamard_soft(i) = mean(BER_hadamard_soft);

end


figure(1);
semilogy(snr,BER_SNR);
hold on;
semilogy(snr+log10(code_rate),BER_SNR_hadamard_hard);
semilogy(snr+log10(code_rate),BER_SNR_hadamard_soft);
legend('Without Channel Coded','Hadamard Hard Decoded','Hadamard Soft Decoded');
grid on;
xlabel('SNR in dB')
ylabel('Bit Error Rate (BER)')


















