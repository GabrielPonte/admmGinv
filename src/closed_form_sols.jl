function closed_form1n(Y,rho_inv)
    return sign.(Y).*max.(abs.(Y) .- (rho_inv), 0)
end

function closed_form21n(Y,rho_inv)
    # Initialize E with zeros
    E = zeros(eltype(Y), size(Y))
    # Compute norms of each row of Y
    norm_rows = [norm(Y[i, :]) for i in axes(Y,1)]
    # Compute norms of each row of Y directly within the condition
    row_idx = rho_inv .< norm_rows
    # Compute the multiplication factor for selected rows
    mult_factor = (norm_rows[row_idx] .- rho_inv) ./ norm_rows[row_idx]
    # Multiply the factors of the selected rows of Y and assign to E
    E[row_idx,:] =  mult_factor .* Y[row_idx,:]
    return E
end

function closed_form2120n(Y,rho_inv,n20)
    # Initialize E with zeros
    E = zeros(eltype(Y), size(Y))
    # Compute norms of each row of Y
    norm_rows = [norm(Y[i, :]) for i in axes(Y,1)]
    # Compute norms of each row of Y directly within the condition
    row_idx0 = partialsortperm(norm_rows, 1:n20, rev=true);
    last_idx = findlast(norm_rows[row_idx0] .> rho_inv);
    if isnothing(last_idx)
        row_idx = [];
    else
        row_idx  = row_idx0[(1:last_idx)];
    end
    # Compute the multiplication factor for selected rows
    mult_factor = (norm_rows[row_idx] .- rho_inv) ./ norm_rows[row_idx]
    # Multiply the factors of the selected rows of Y and assign to E
    E[row_idx,:] =  mult_factor .* Y[row_idx,:]
    return E
end

function closed_form20n(Y, n20)
    # Initialize E with zeros
    E = zeros(eltype(Y), size(Y))
    # Compute norms of each row of Y
    norm_rows = [norm(Y[i, :]) for i in axes(Y,1)]
    # Get indices of top n20 rows based on norm_rows
    row_idx = partialsortperm(norm_rows, 1:n20, rev=true)
    # Assign rows of Y to E based on row_idx
    @views E[row_idx, :] .= Y[row_idx, :]
    return E
end