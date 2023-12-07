function sta = calculateGratingFlickerSTA(stimmat, runningspikes, stimorder)


[Nt, Nstimuli] = size(stimorder);
Ncells         = size(runningspikes, 1);
[Nyx, Nbase]   = size(stimmat);


% use GPU in a smart way: break in chunks of different pixels!!!

%--------------------------------------------------------------------------
%calculate stim chunks
g=gpuDevice(1); availableMem = g.AvailableMemory*8 - 5e9; %in bits with buffer
staMem    = Ncells*Nyx*32; 
spikesMem = Ncells*(Nstimuli)*32; 
stimMem   = Nyx*Nstimuli*32;
totalMem  = (staMem+spikesMem+stimMem);
Nchunks   = ceil(totalMem/availableMem); 
chunkSize = floor(Nyx/Nchunks);
%--------------------------------------------------------------------------
sta      = zeros(Ncells,Nyx,Nt,'single'); %preallocate sta in RAM
spikes   = single(gpuArray(runningspikes));
%--------------------------------------------------------------------------

msg = [];tic;
for ichunk = 1:Nchunks
    dimstart  = (ichunk - 1) * chunkSize + 1;
    dimstop   = min(ichunk * chunkSize, Nyx);
    chunkdims = dimstop - dimstart + 1;
    gpusta    = gpuArray.zeros(Ncells, chunkdims, Nt,'single'); %preallocate sta in RAM
    chunkstim = single(gpuArray(stimmat(dimstart:dimstop, :)));
    for it=1:Nt
        gpusta(:,:,it) = gpusta(:,:,it) + spikes * chunkstim(:, stimorder(it, :))';
    end 
    sta(:, dimstart:dimstop, :) = gather(gpusta);
    %--------------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Dimension chunk %d/%d. Time elapsed %2.2f s...\n', ...
        ichunk, Nchunks,toc);
    fprintf(msg);
    %--------------------------------------------------------------------------
    clear gpusta;

end
sumspikes = sum(runningspikes, 2);
sta       = bsxfun(@rdivide, sta, sumspikes);%cast sta to double and divide by nspikes


% 
% 

% for istim = 1:Nstim
%     currorder = stimorder{istim};
%     spikesbin = single(runningspikes{istim}(:, Nt:end, :));
%     Nblocks   = size(currorder, 2);
%     
%     for iblock = 1:Nblocks
%         spikes   = spikesbin(:, :, iblock);
%         stimulus = stimmat(:, currorder(:, iblock));
%         stimulus = 2*single(stimulus)/255-1;%making stimulus into single mat of 1s and -1s
%         
%         for it=1:Nt
%             sta(:,:,it) = sta(:,:,it) + spikes * stimulus(:, it:(end-Nt+it))';
%         end
%         %--------------------------------------------------------------------------
%         fprintf(repmat('\b', 1, numel(msg)));
%         msg = sprintf('Stim %d/%d, block %d/%d. Time elapsed %2.2f s...\n', ...
%             istim, Nstim, iblock, Nblocks,toc);
%         fprintf(msg);
%         %--------------------------------------------------------------------------
%     end
%     
%     
% end
% 
% 
% 
% 
% %--------------------------------------------------------------------------
% 
% 
% spikesbin = single(spikesbin(:,Nt:frameN,:)); %casting spikes to single
% sumspikes = sum(spikesbin(:,:), 2);
% %--------------------------------------------------------------------------
% %calculate cell chunks
% g=gpuDevice(1); availableMem = g.AvailableMemory*8 - 8e9; %in bits with buffer
% staMem    = cellN*Nyx*Nt*32; 
% spikesMem = cellN*(frameN-Nt+1)*32; 
% stimMem   = Nyx*frameN*32;
% totalMem  = (staMem+spikesMem+stimMem); %buffer of 500MB added
% Nchunks   = ceil(totalMem/availableMem); chunkSize=floor(cellN/Nchunks);
% %--------------------------------------------------------------------------
% sta      = zeros(cellN,Nyx,Nt,'single'); %preallocate sta in RAM
% stimulus = zeros(Nyx,frameN,'single','gpuArray'); %preallocate stimulus
% %--------------------------------------------------------------------------
% msg = []; tic;
% for iChunk=1:Nchunks
%     cellstart=(iChunk-1)*chunkSize+1;
%     cellchunk=min(chunkSize, cellN-cellstart+1);
%     cellend=cellstart+cellchunk-1;
%     chunkspikesbin=spikesbin(cellstart:cellend,:,:); %casting spikes to single
% 
%     chunksta = zeros(cellchunk, Nyx,    Nt, 'single','gpuArray'); %preallocate sta
%     spikes   = zeros(cellchunk, frameN-Nt+1, 'single','gpuArray');
%     seed     = oseed;
%     %--------------------------------------------------------------------------
%     for iBlock=1:blockN 
%         [stimulus(:),seed]=ran1(seed,Nyx*frameN);
%         stimulus  = 2*single(stimulus>0.5)-1;%making stimulus into single mat of 1s and -1s
%         spikes(:) = chunkspikesbin(:,:,iBlock);
%         for it=1:Nt
%             chunksta(:,:,it) = chunksta(:,:,it) + ...
%                 contrast * spikes * stimulus(:, it:(end-Nt+it))';
%         end
%     end
%     %--------------------------------------------------------------------------
%     sta(cellstart:cellend,:,:)=gather(chunksta);
%     clear chunksta;
%     %--------------------------------------------------------------------------
%     fprintf(repmat('\b', 1, numel(msg)));
%     msg = sprintf('Chunk %d/%d. Time elapsed %2.2f s...\n', iChunk,Nchunks,toc);
%     fprintf(msg);
% end
% sta = bsxfun(@rdivide, sta, sumspikes);%cast sta to double and divide by nspikes
end