function isListening = is_listening(vad, channel)
    % Validate the channel number
    numChannels = size(vad,2);
    if channel < 1 || channel > numChannels
        error('Channel number is out of bounds.');
    end
    
    % Invert the voice activity detection for the specified channel
    noActivityInChannel = ~vad(:, channel);
    
    % Get all channels except the specified one
    allChannels = 1:numChannels;
    otherChannels = allChannels ~= channel;
    
    % Check if there's any voice activity in the other channels
    anyActivityInOtherChannels = any(vad(:, otherChannels), 2);
    
    % Determine if the system is listening
    isListening = noActivityInChannel & anyActivityInOtherChannels;
end