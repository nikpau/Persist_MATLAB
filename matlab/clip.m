function x = clip(x,lowerBound,upperBound)

    x = max(min(upperBound,x),lowerBound);
    
end

