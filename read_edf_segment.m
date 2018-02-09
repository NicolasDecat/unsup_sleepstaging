function datamat = read_edf_segment(dataset,fs,interval,nChannel)
% ================================================
% Segment data into certain time interval after reading data from edf file
% 
% Input: dataset--matrix of data to be segmented 
%        fs--sampling frequency (Hz)
%        totaltime--total time of the data set (seconds)
%        interval--time duration for each segment (seconds)
%        nChannel-- number of channels in the dataset
% Output: data_mat--matrix of data in interval
%
% ===================================================
    %% Check input argument
    if nargin < 1
        error('Insufficient input\n');
    end

    %% Determine length of segment
    [expChannel,nSample] = size(dataset);

    % Check if the expected number of channel matched with specified channel
    if (expChannel~=nChannel)
        error('Number of channel specified does not match number of channel in dataset.\n')
    end

    % Length of each segment
    LSeg = interval*fs;

    % Number of segment
    nSeg = ceil(nSample/LSeg);

    % Size of final data matrix
    n_row = nChannel*nSeg;
    n_column = LSeg;
    datamat = ones(n_row,n_column)*NaN;

    %% Segment into interval (Lseg samples/interval)
    tempmat = zeros(nSeg,LSeg);     % Temporary matrix

    for k=1:nChannel  
        n=1;
        for i=1:nSeg
            tempmat(i,:) = dataset(k,n:n+(LSeg-1));
            n=n+LSeg;
        end
        if (k==1)
            datamat = tempmat;
        else
            datamat = [datamat;tempmat];
        end
    end 
end