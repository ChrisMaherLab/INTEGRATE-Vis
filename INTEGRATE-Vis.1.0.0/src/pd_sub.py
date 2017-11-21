import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import gene_list 
from scipy.special import stdtr

write_sort = True

# plot set up
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
sample_list = []
oncogene_max_exp = 0
partner_max_exp = 0
x_coor = 1.2
offset = 1.5
bar_height = 1
bar_width = 0.1
font_size = 10.5

#vertial parameter for different size
y_factor = 1.0

def p_value(lista, listb):
    na=len(lista)
    nb=len(listb)
    adof = na - 1
    bdof = nb - 1
    if adof == 0:
        lista=lista+lista
        adof=1
        na=2
    if bdof == 0:
        listb=listb+listb
        bdof=1
        nb=2
    abar=np.mean(lista)
    bbar=np.mean(listb)
    avar = np.var(lista)
    bvar = np.var(listb)
    tf = (abar - bbar) / np.sqrt(avar/na + bvar/nb)
    dof = (avar/na + bvar/nb)**2 / (avar**2/(na**2*adof) + bvar**2/(nb**2*bdof))
    pf = 2*stdtr(dof, -np.abs(tf))
    return pf

def get_suffix(name):
    if name in gene_list.oncogene_list:
        return " (oncogene)"
    if name in gene_list.tumor_suppressor_list:
        return " (tumor suppressor)" 
    return ""

class sample:
    name = ""
    three_prime_exp = 0
    five_prime_exp = 0
    tumor = True
    fusion = False
    exp_text = "Expression"
    y_mode = 'independent'
    max_value = -10000.0
    # string var annotating both genes

    def __init__(self, name, tumor, fusion, three_prime_exp, five_prime_exp, exp_text, y_mode, max_value):
        self.name = name
        self.five_prime_exp = five_prime_exp
        self.three_prime_exp = three_prime_exp
        self.tumor = tumor
        self.fusion = fusion
        self.exp_text = exp_text
        self.y_mode = y_mode
        self.max_value = float(max_value)
        self.write_sort=write_sort        

    def __repr__(self):
        return str(self.name) + " " + str(self.tumor) + " " + str(self.fusion) + " " + str(self.three_prime_exp) + " " + str(self.five_prime_exp) + "\n"

    def draw(self, x_coor):
        if (self.tumor and self.fusion):
            color = "#252525"
        elif (self.tumor):
            color = "#ae017e"
        else:
            color = "#006d2c"
# here bar_width * something; 1.3 is the global space for separating bars
        maxa=1
        maxb=1
        if self.y_mode != "independent":
            maxa=max(oncogene_max_exp,partner_max_exp)
            maxb=max(oncogene_max_exp,partner_max_exp)
        else:
            maxa=oncogene_max_exp
            maxb=partner_max_exp
        value_a = sample.three_prime_exp
        if value_a>maxa:
            value_a=maxa
        value_b = sample.five_prime_exp
        if value_b>maxb:
            value_b = maxb

        oncogene_bar = patches.Rectangle(
            (x_coor, 0.5 + 1.3 * bar_width), bar_width, bar_height * (value_a/maxa),
            facecolor= "#cb181d", ec = "none"
        )
        oncogene_base = patches.Rectangle(
            (x_coor - bar_width*0.15, 0.5 - bar_width*y_factor) , bar_width*1.3, bar_width*y_factor,
            facecolor= color, ec = "none"
        )
        partner_bar = patches.Rectangle(
            (x_coor, 0.5 + 1.3 * bar_width + offset), bar_width, bar_height * (value_b/maxb),
            facecolor= "blue", ec = "none"
        )
        partner_base = patches.Rectangle(
            (x_coor - bar_width*0.15, 0.5 + offset - bar_width*y_factor), bar_width*1.3, bar_width*y_factor,
            facecolor= color, ec = "none"
        )

        ax.add_patch(oncogene_bar)
        ax.add_patch(oncogene_base)
        ax.add_patch(partner_bar)
        ax.add_patch(partner_base)

# use three lists

def three_prime_sort_key(s):
    if (not s.tumor):
        return 0 + s.three_prime_exp
    elif (s.fusion):
        return 2*(10**9) + s.three_prime_exp
    else:
        return (10**9) + s.three_prime_exp

def five_prime_sort_key(s):
    if (not s.tumor):
        return 0 + s.five_prime_exp
    elif (s.fusion):
        return 2*(10**9) + s.five_prime_exp
    else:
        return (10**9) + s.five_prime_exp

def default_sort(sample_list):#here onco means three prime, partner means 5 prime
    oncogene_exp_normal = []
    oncogene_exp_tumor = []
    oncogene_exp_fusion = []
    partner_exp_normal = []
    partner_exp_tumor = []
    partner_exp_fusion = []
    for sample in sample_list:
        if(not sample.tumor):
            oncogene_exp_normal.append(sample.three_prime_exp)
            partner_exp_normal.append(sample.five_prime_exp)
        elif(sample.fusion):
            oncogene_exp_fusion.append(sample.three_prime_exp)
            partner_exp_fusion.append(sample.five_prime_exp)
        else:
            oncogene_exp_tumor.append(sample.three_prime_exp)
            partner_exp_tumor.append(sample.five_prime_exp)

    if (len(oncogene_exp_normal) != 0):
        pval3=p_value(oncogene_exp_tumor+oncogene_exp_fusion,oncogene_exp_normal)
        pval5=p_value(partner_exp_tumor+partner_exp_fusion,partner_exp_normal)
        if pval3<pval5:
#       if (abs((np.mean(oncogene_exp_tumor)+1.0)/ (np.mean(oncogene_exp_normal)+1.0)) > abs((np.mean(partner_exp_tumor)+1.0) / (np.mean(partner_exp_normal)+1.0))):
            sample_list.sort(key = three_prime_sort_key)
            return sample_list #are arguments passed by reference?? If so can simplify
        else:
            sample_list.sort(key = five_prime_sort_key)
            return sample_list
    else:
        #print np.median(oncogene_exp_fusion),np.mean(oncogene_exp_fusion)
        #print np.median(oncogene_exp_tumor),np.mean(oncogene_exp_tumor)
        #print np.median(partner_exp_fusion),np.mean(partner_exp_fusion)
        #print np.median(partner_exp_tumor),np.mean(partner_exp_tumor)
        pval3=p_value(oncogene_exp_tumor,oncogene_exp_fusion)
        pval5=p_value(partner_exp_tumor,partner_exp_fusion)
        #print "pval3,5=",pval3,pval5

        #if (abs((np.mean(oncogene_exp_fusion)+1.0) / (np.mean(oncogene_exp_tumor)+1.0)) > abs((np.mean(partner_exp_fusion)+1.0) / (np.mean(partner_exp_tumor))+1.0)):
        if pval3 < pval5:
            sample_list.sort(key = three_prime_sort_key)
            return sample_list
        else:
            sample_list.sort(key = five_prime_sort_key)
            return sample_list


# should be given as argument
def visualize(pd_matrix_dir, five_prime_name, three_prime_name, sort_mode, output_dir, exp_text, y_mode, max_value, write_sort):
    global sample, sample_list, oncogene_max_exp
    global partner_max_exp
    global offset
    global x_coor
    global ax
    global fig
    global bar_height
    global bar_width
    global font_size
    # read in data
    with open(pd_matrix_dir) as input_file: # add more binary cols indicating kinase (y/n) oncogene (y/n) ...
        next(input_file)
        for line in input_file:
            #print(line)
            tmp = line.split("\t")
            name = tmp[0]
            tumor = bool(int(tmp[1]))
            fusion = bool(int(tmp[2]))
            three_prime_exp = float(tmp[3])
            five_prime_exp = float(tmp[4])
            sample_list.append(sample(name, tumor, fusion, three_prime_exp, five_prime_exp, exp_text, y_mode, max_value))

    #calculate y_factor
    tt=len(sample_list)+0.0
    global y_factor
    y_factor=tt*tt/100000.0+tt/50.0+1
    #print "###########  y_factor=",y_factor

    # sort the expression
    if (sort_mode == 0): # default sorting
         sample_list = default_sort(sample_list)
    elif (sort_mode == 5):
        sample_list.sort(key = five_prime_sort_key)
    elif (sort_mode == 3):
        sample_list.sort(key = three_prime_sort_key)

    # print
    if write_sort:
        print ("%s\t%s\t%s\t%s\t%s" % ("sample name","is tumor", "is fusion", "5p exp", "3p exp"))
        for x in range(len(sample_list)):
            print ("%s\t%s\t%s\t%s\t%s" % (sample_list[x].name, str(sample_list[x].tumor), str(sample_list[x].fusion), str(sample_list[x].five_prime_exp), str(sample_list[x].three_prime_exp)))

    bar_height = 1.5
    bar_width = 4.0/len(sample_list)
    #print("bar_width")
    #print(bar_width)
    oncogene_max_exp = max(sample.three_prime_exp for sample in sample_list)
    partner_max_exp = max(sample.five_prime_exp for sample in sample_list)
    offset = bar_height * 1.5

    #adjust to max_value
    
    if float(max_value) > 0.0:
        oncogene_max_exp=float(max_value)
        partner_max_exp=float(max_value)

    # draw the expression bars
    for sample in sample_list:
        sample.draw(x_coor)
        x_coor += bar_width * 1.3;

    xlim = x_coor
    ax.set_xlim([0, xlim])
    ax.set_ylim([0, 2.1 * bar_height + offset])

    # draw tick
    maxa=1
    maxb=1
    if y_mode != "independent":
        maxa=max(oncogene_max_exp,partner_max_exp)
        maxb=max(oncogene_max_exp,partner_max_exp)
    else:
        maxa=oncogene_max_exp
        maxb=partner_max_exp

    plt.text(1.0,  0.5+bar_height*1.0, str(maxa), fontsize = font_size*0.75, verticalalignment ='bottom', horizontalalignment='right', color = 'red')
    plt.text(1.1,  0.5+bar_height*1.0, "_", fontsize = font_size, verticalalignment ='bottom', color = 'red')
    plt.text(1.0,  0.5, "0", fontsize = font_size*0.75, verticalalignment ='bottom', horizontalalignment='right', color = 'red')
    plt.text(1.1,  0.5, "_", fontsize = font_size, verticalalignment ='bottom',color = 'red')
   
    plt.text(1.0,  0.5+bar_height*1.0 + offset, str(maxb), fontsize = font_size*0.75, verticalalignment ='bottom', horizontalalignment='right', color = 'red')
    plt.text(1.1,  0.5+bar_height*1.0 + offset, "_", fontsize = font_size, verticalalignment ='bottom',color = 'red')
    plt.text(1.0,  0.5+offset, "0", fontsize = font_size*0.75, verticalalignment ='bottom', horizontalalignment='right',color = 'red')
    plt.text(1.1,  0.5+offset, "_", fontsize = font_size, verticalalignment ='bottom',color = 'red')

    # legends
    plt.text(-0.2, 1.0 - bar_width*y_factor, exp_text, fontsize = font_size, horizontalalignment = 'left')
    plt.text(-0.2, 1.0 - bar_width*y_factor + offset, exp_text, fontsize = font_size, horizontalalignment = 'left')
    plt.text(-0.2, 0.5 - bar_width*y_factor, "Type", fontsize = font_size, horizontalalignment = 'left')
    plt.text(-0.2, 0.5 - bar_width*y_factor + offset, "Type", fontsize = font_size, horizontalalignment = 'left')
    plt.text(xlim/2+0.6, 0.625 + bar_height * 1.1, three_prime_name + get_suffix(three_prime_name), horizontalalignment = 'center', fontsize = font_size)
    plt.text(xlim/2+0.6, 0.625 + bar_height * 1.1 + offset, five_prime_name+get_suffix(five_prime_name), horizontalalignment = 'center', fontsize = font_size)

    #y_factor*1.0 means on x axis
    legend_x_coor = xlim/2 + 0.8 - 1.3
    plt.text(legend_x_coor - 1 + bar_width*y_factor*1.0, 0, "Normal", fontsize = font_size)
    #legend need to change
    color_legend1 = patches.Rectangle((legend_x_coor - 1, 0), bar_width*y_factor*0.5, bar_width*y_factor*1.2, facecolor= "#006d2c", ec = "none")
    plt.text(legend_x_coor + bar_width*y_factor*1.0, 0, "Tumor without Fusion", fontsize = font_size)
    color_legend2 = patches.Rectangle((legend_x_coor, 0), bar_width*y_factor*0.5, bar_width*y_factor*1.2, facecolor= "#ae017e", ec = "none")
    plt.text(legend_x_coor + 2 + bar_width*y_factor*1.0, 0, "Tumor with Fusion", fontsize = font_size)
    color_legend3 = patches.Rectangle((legend_x_coor + 2, 0), bar_width*y_factor*0.5, bar_width*y_factor*1.2, facecolor= "#252525", ec = "none")
    color_legends = [color_legend1, color_legend2, color_legend3]
    for p in color_legends:
        ax.add_patch(p)
    plt.axis('off')
    #print(output_dir)
    fig.tight_layout()
    plt.savefig(output_dir,bbox_inches='tight')

