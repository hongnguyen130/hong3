function create_channel_dataset_SR(Ns, parGenCh, save_name)
    % Initialize the channel dataset
    rng(2025);%de tai lap, co the doi seed
    samples=repmat(struct(),Ns,1);
    for i = 1:Ns
        fprintf('Generating sample %d/%d\n',i,Ns);
        % ----- Generate small-scale fading -----
        ch = generate_channel_SR(parGenCh);
        % ----- Save per-sample data -----
        samples(i).hBI = ch.hBI;
        samples(i).g1  = ch.g1;
        samples(i).g2  = ch.g2;
        samples(i).gE  = ch.gE;

        samples(i).hB1 = ch.hB1;
        samples(i).hB2 = ch.hB2;
        samples(i).hBE = ch.hBE;
        
    end
    channel_dataset.Ns = Ns;
    channel_dataset.samples = samples;
    save(save_name, 'channel_dataset', '-v7.3');
    fprintf('Saved dataset to %s\n', save_name);
end