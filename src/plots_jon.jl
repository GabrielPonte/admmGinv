using Plots, CSV,DataFrames,LaTeXStrings


# Assuming your CSV file is named "data.csv"
data = CSV.read(".//csv//1norm.csv",DataFrame,header=false, normalizenames=false)

# Assuming your CSV file has columns named "x" and "y"
x_vals = [100,200,300,400,500,1000,2000,3000,4000,5000]
# x = data[!, 3]
y = data[!, 3]
y_vals = [data[!,1] data[!,2] data[!,3] data[!,4]]
str1 = L"\textrm{ADMM}_{1}^\epsilon"
str2 = L"\textrm{ADMM}_{1}^{}"
my_legend = [L"{P}^{1}_{123}" L"\mathcal{P}^{1^{^{}}}_{123}" str1 str2]
# latexstring("\$\\textrm{LS}\$ \n \$\\textrm{[ $(time[8]) ]}\$")
plot(
        x_vals,
        y_vals,
        lw = 1.25,
        yaxis=:log,
        legend=:bottomright,
        marker=:circle,
        xticks =  (000:1000:5000),
        yticks =  10.0.^(5:-1:-15),
        xlabel = L"m",
        ylabel = "Elapsed time (sec)",
        label = my_legend,
        legendfontsize=11.0,
        # legendfonthalign = :right,
        palette=[palette(:darkrainbow)[5];palette(:default)[[1,3,4]]],
        # palette = palette(:darkrainbow)#[[1,3]]
    )

    # savefig("plot_1norm_v4.pdf")


###############################
# Assuming your CSV file is named "data.csv"
data = CSV.read(".//csv//21norm.csv",DataFrame,header=false, normalizenames=false)
str1 = L"\textrm{ADMM}_{2,1}^\epsilon"
str2 = L"\textrm{ADMM}_{2,1}^{}"
my_legend = [L"P^{2,1}_1" L"\mathcal{P}^{2,1^{^{}}}_{1}" str1 str2]

# data = CSV.read(".//csv//1norm.csv",DataFrame,header=false, normalizenames=false)
# str1 = L"\textrm{ADMM}_{1}^\epsilon"
# str2 = L"\textrm{ADMM}_{1}^{}"
# my_legend = [L"P_{1}" L"\mathcal{P}_{1}" str1 str2]

x_vals = [100,200,300,400,500,1000,2000,3000,4000,5000]
y_vals = [data[!,1] data[!,2] data[!,3] data[!,4]]


    plot(
        x_vals,
        y_vals,
        lw = 1.25,
        yaxis=:log,
        legend=:bottomright,
        marker=:circle,
        xticks =  (000:1000:5000),
        yticks =  10.0.^(5:-1:-15),
        xlabel = L"m",
        ylabel = "Elapsed time (sec)",
        label = my_legend,
        legendfontsize=11.0,
        # legendfonthalign = :right,
        palette=[palette(:darkrainbow)[5];palette(:default)[[1,3,4]]],
        # palette = palette(:darkrainbow)#[[1,3]]
    )
    # savefig("plot_21norm_v4.pdf")


##############################################################################

# PLOT 2,0  - L4
data = CSV.read(".//csv//20_L4.csv",DataFrame,header=false, normalizenames=false)
x_vals = [100,200,300,400,500,1000,2000,3000,4000,5000]
y_vals = [0.5*data[!,1] data[!,2]]
y_vals45 = [0.5*data[!,4] data[!,5]]
time = data[!,3]
str1 = L"0.5\;{\|\|H\|\|}_{2,1}^{}"
str2 = L"{\|\|H\|\|}_{2,0}"
my_legend = [str1 str2]


str_x1 = latexstring("\$\\textrm{ADMM}_{2,1}\$ \n \$\\textrm{[ $(time[1]) ]}\$")
str_x2 = latexstring("\$0.25\$ \n \$\\textrm{[ $(time[2]) ]~\\;\\,}\$")
str_x3 = latexstring("\$0.50\$ \n \$\\textrm{[ $(time[3]) ]~\\;\\,}\$")
str_x4 = latexstring("\$0.75\$ \n \$\\textrm{[ $(time[4]) ]~\\;\\,}\$")
str_x5 = latexstring("\$0.80\$ \n \$\\textrm{[ $(time[5]) ]~\\;\\,}\$")
str_x6 = latexstring("\$0.90\$ \n \$\\textrm{[ $(time[6]) ]~\\;\\,}\$")
str_x7 = latexstring("\$0.95\$ \n \$\\textrm{[ $(time[7]) ]~\\;\\,}\$")
str_x8 = latexstring("\$\\textrm{LS}\$ \n \$\\textrm{[ $(time[8]) ]~\\;\\,}\$")

x_vals = [str_x1, str_x2, str_x3,str_x4,str_x5,str_x6,str_x7, str_x8]

scatter(x_vals,y_vals45,palette = palette(:default)[[1,2]],label = nothing)
plot!(
    x_vals,
    y_vals,
    lw = 1.25,
    label = my_legend,
    xlabel = L"\omega",
    ylabel = L"\|\|H\|\|",
    legendfontsize=11.0,
    xtickfontsize=10,
    legend=:top
)
savefig("plot_20norm_L4_v6.pdf")


########################################
# PLOT 2,0  - L5
data = CSV.read(".//csv//20_L5.csv",DataFrame,header=false, normalizenames=false)
x_vals = [100,200,300,400,500,1000,2000,3000,4000,5000]
data1 = copy(data[!,1]); 
data2 = copy(data[!,2]); 
y_vals = [0.5*data[!,1] data[!,2]]
y_vals45 = [0.5*data[!,4] data[!,5]]
time = data[!,3]
str1 = L"0.5\;{\|\|H\|\|}_{2,1}^{}"
str2 = L"{\|\|H\|\|}_{2,0}"
my_legend = [str1 str2]

str_x1 = latexstring("\$\\textrm{ADMM}_{2,1}\$ \n \$\\textrm{[ $(time[1]) ]}\$")
str_x2 = latexstring("\$0.25\$ \n \$\\textrm{[ $(time[2]) ]~\\;\\,}\$")
str_x3 = latexstring("\$0.50\$ \n \$\\textrm{[ $(time[3]) ]~\\;\\,}\$")
str_x4 = latexstring("\$0.75\$ \n \$\\textrm{[ $(time[4]) ]~\\;\\,}\$")
str_x5 = latexstring("\$0.80\$ \n \$\\textrm{[ $(time[5]) ]~\\;\\,}\$")
str_x6 = latexstring("\$0.90\$ \n \$\\textrm{[ $(time[6]) ]~\\;\\,}\$")
str_x7 = latexstring("\$0.95\$ \n \$\\textrm{[ $(time[7]) ]~\\;\\,}\$")
str_x8 = latexstring("\$\\textrm{LS}\$ \n \$\\textrm{[ $(time[8]) ]~\\;\\,}\$")

x_vals = [str_x1, str_x2, str_x3,str_x4,str_x5,str_x6,str_x7, str_x8]


scatter(x_vals,y_vals45,palette = palette(:default)[[1,2]],label = nothing)
plot!(
    x_vals,
    y_vals,
    lw = 1.25,
    label = my_legend,
    xlabel = L"\omega",
    ylabel = L"\|\|H\|\|",
    legendfontsize=11.0,
    xtickfontsize=10,
    legend=:top,
    yticks =  (950:350:2400),
)
savefig("plot_20norm_L5_v6.pdf")