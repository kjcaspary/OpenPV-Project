# -*- coding: utf-8 -*-

#from __future__ import (absolute_import, division)

import pandas    
import numpy

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.colors import rgb2hex, Normalize
from matplotlib.patches import Polygon
from matplotlib.colorbar import ColorbarBase

from bokeh.plotting import figure, show, output_file
from bokeh.palettes import Magma256 as palette
from bokeh.models import (
    ColumnDataSource,
    HoverTool,
    LogColorMapper
)


def StateConvert(states_in):
    states_out = states_in[:]
    states = {
            'AK': 'Alaska',
            'AL': 'Alabama',
            'AR': 'Arkansas',
            'AS': 'American Samoa',
            'AZ': 'Arizona',
            'CA': 'California',
            'CO': 'Colorado',
            'CT': 'Connecticut',
            'DC': 'District of Columbia',
            'DE': 'Delaware',
            'FL': 'Florida',
            'GA': 'Georgia',
            'GU': 'Guam',
            'HI': 'Hawaii',
            'IA': 'Iowa',
            'ID': 'Idaho',
            'IL': 'Illinois',
            'IN': 'Indiana',
            'KS': 'Kansas',
            'KY': 'Kentucky',
            'LA': 'Louisiana',
            'MA': 'Massachusetts',
            'MD': 'Maryland',
            'ME': 'Maine',
            'MI': 'Michigan',
            'MN': 'Minnesota',
            'MO': 'Missouri',
            'MP': 'Northern Mariana Islands',
            'MS': 'Mississippi',
            'MT': 'Montana',
            'NA': 'National',
            'NC': 'North Carolina',
            'ND': 'North Dakota',
            'NE': 'Nebraska',
            'NH': 'New Hampshire',
            'NJ': 'New Jersey',
            'NM': 'New Mexico',
            'NV': 'Nevada',
            'NY': 'New York',
            'OH': 'Ohio',
            'OK': 'Oklahoma',
            'OR': 'Oregon',
            'PA': 'Pennsylvania',
            'PR': 'Puerto Rico',
            'RI': 'Rhode Island',
            'SC': 'South Carolina',
            'SD': 'South Dakota',
            'TN': 'Tennessee',
            'TX': 'Texas',
            'UT': 'Utah',
            'VA': 'Virginia',
            'VI': 'Virgin Islands',
            'VT': 'Vermont',
            'WA': 'Washington',
            'WI': 'Wisconsin',
            'WV': 'West Virginia',
            'WY': 'Wyoming'
            }
#    print states_in,states_in[0:2]
#    print states[states_in[0]]
    for i in range(len(states_in)):
        
        states_out[i] = states[states_in[i]]
#        print states_in[i],states_out[i]
    return states_out
   
def BokehMap(dict_data,d_min, d_max,title_str,cb_str):
    us_states = numpy.load('us_states.npy').item()
    
#    us_counties = us_counties.data.copy()
#    unemployment = unemployment.data
    
    del us_states["HI"]
    del us_states["AK"]
    #del dict_data["HI"]
    #del dict_data["AK"]
    
    state_xs = [us_states[code]["lons"] for code in us_states]
    state_ys = [us_states[code]["lats"] for code in us_states]
    
    #county_xs=[us_counties[code]["lons"] for code in us_counties if us_counties[code]["state"] not in ["ak", "hi", "pr", "gu", "vi", "mp", "as"]]
    #county_ys=[us_counties[code]["lats"] for code in us_counties if us_counties[code]["state"] not in ["ak", "hi", "pr", "gu", "vi", "mp", "as"]]
    
    #colors = ["#F1EEF6", "#D4B9DA", "#C994C7", "#DF65B0", "#DD1C77", "#980043"]
    
    #county_colors = []
    state_colors = []
    state_data = []
    for state_id in us_states.keys():
    #    state_id = state_id.lower()
    #    if us_counties[county_id]["state"] in ["ak", "hi", "pr", "gu", "vi", "mp", "as"]:
    #        continue
        try:
    #        rate = numpy.nanmean([unemployment[code] for code in us_counties
    #                             if us_counties[code]['state']==state_id])
            data = dict_data[state_id]
            print state_id, data
            state_data.append(data)
    #        state_colors.append(colors[idx])
        except KeyError:
    #        state_colors.append("black")
            state_data.append(0.0)
    
    output_file("choropleth.html", title="choropleth.py example")
    
    #palette.reverse()
    color_mapper = LogColorMapper(palette=palette)
    
    source = ColumnDataSource(data=dict(
        x=state_xs,
        y=state_ys,
        name=us_states.keys(),
        data=state_data,
    ))
    #p = figure(title="US Unemployment 2009", toolbar_location="left",
    #    plot_width=1100, plot_height=700)
    TOOLS = "pan,wheel_zoom,reset,hover,save"


    p = figure(
        title=title_str, tools=TOOLS,
        plot_width=1100, plot_height=700
        )
    #p.patches(county_xs, county_ys, fill_color=county_colors, fill_alpha=0.7,
    #    line_color="white", line_width=0.5)
    #p.patches(state_xs, state_ys,source=source, fill_color=state_colors,fill_alpha=0.7,
    #    line_color="#884444", line_width=2)
    p.patches('x', 'y', source=source,
              fill_color={'field': 'data', 'transform': color_mapper},
              fill_alpha=0.7, line_color="#884444", line_width=1.0)
              
    hover = p.select_one(HoverTool)
    hover.point_policy = "follow_mouse"
    hover.tooltips = [
        ("Name", "@name"),
        (title_str+' ('+cb_str+')', "@data"),
        ("(Long, Lat)", "($x, $y)"),
        ]
    show(p) 
## adapted from https://github.com/matplotlib/basemap/blob/master/examples/fillstates.py
def PlotMap(dict_data,d_min, d_max,title_str,cb_str):

    
    fig, ax = plt.subplots()
    
    # Lambert Conformal map of lower 48 states.
    m = Basemap(llcrnrlon=-119,llcrnrlat=20,urcrnrlon=-64,urcrnrlat=49,
                projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
    
    # Mercator projection, for Alaska and Hawaii
    m_ = Basemap(llcrnrlon=-190,llcrnrlat=20,urcrnrlon=-143,urcrnrlat=46,
            projection='merc',lat_ts=20)  # do not change these numbers

    ## http://www.census.gov/geo/www/cob/st2000.html
    shp_info = m.readshapefile('st99_d00','states',drawbounds=True,
                               linewidth=0.45,color='gray')
    shp_info_ = m_.readshapefile('st99_d00','states',drawbounds=False)
    
    ## population density by state from
    ## http://en.wikipedia.org/wiki/List_of_U.S._states_by_population_density
    
    plot_data = dict_data

    colors={}
    statenames=[]
    cmap = plt.cm.hot_r # use 'reversed hot' colormap
    vmin = d_min; vmax = d_max # set range.
    norm = Normalize(vmin=vmin, vmax=vmax)
    for shapedict in m.states_info:
        statename = shapedict['NAME']
        # skip DC and Puerto Rico.
        if statename not in ['District of Columbia','Puerto Rico']:
            dat = plot_data[statename]
            # calling colormap with value between 0 and 1 returns
            # rgba value.  Invert color range (hot colors are high
            # population), take sqrt root to spread out colors more.
            colors[statename] = cmap(np.sqrt((dat-vmin)/(vmax-vmin)))[:3]
        statenames.append(statename)

    for nshape,seg in enumerate(m.states):
        # skip DC and Puerto Rico.
        if statenames[nshape] not in ['Puerto Rico', 'District of Columbia']:
            color = rgb2hex(colors[statenames[nshape]])
            poly = Polygon(seg,facecolor=color,edgecolor=color)
            ax.add_patch(poly)
        
    AREA_1 = 0.005  # exclude small Hawaiian islands that are smaller than AREA_1
    AREA_2 = AREA_1 * 30.0  # exclude Alaskan islands that are smaller than AREA_2
    AK_SCALE = 0.19  # scale down Alaska to show as a map inset
    HI_OFFSET_X = -1900000  # X coordinate offset amount to move Hawaii "beneath" Texas
    HI_OFFSET_Y = 250000    # similar to above: Y offset for Hawaii
    AK_OFFSET_X = -250000   # X offset for Alaska (These four values are obtained
    AK_OFFSET_Y = -750000   # via manual trial and error, thus changing them is not recommended.)
    
    for nshape, shapedict in enumerate(m_.states_info):  # plot Alaska and Hawaii as map insets
        if shapedict['NAME'] in ['Alaska', 'Hawaii']:
            seg = m_.states[int(shapedict['SHAPENUM'] - 1)]
            if shapedict['NAME'] == 'Hawaii' and float(shapedict['AREA']) > AREA_1:
                seg = [(x + HI_OFFSET_X, y + HI_OFFSET_Y) for x, y in seg]
                color = rgb2hex(colors[statenames[nshape]])
            elif shapedict['NAME'] == 'Alaska' and float(shapedict['AREA']) > AREA_2:
                seg = [(x*AK_SCALE + AK_OFFSET_X, y*AK_SCALE + AK_OFFSET_Y)\
                       for x, y in seg]
                color = rgb2hex(colors[statenames[nshape]])
            poly = Polygon(seg, facecolor=color, edgecolor='gray', linewidth=.45)
            ax.add_patch(poly)
        
    ax.set_title(title_str+'by State')

    light_gray = [0.8]*3  # define light gray color RGB
    x1,y1 = m_([-190,-183,-180,-180,-175,-171,-171],[29,29,26,26,26,22,20])
    x2,y2 = m_([-180,-180,-177],[26,23,20])  # these numbers are fine-tuned manually
    m_.plot(x1,y1,color=light_gray,linewidth=0.8)  # do not change them drastically
    m_.plot(x2,y2,color=light_gray,linewidth=0.8)

    ax_c = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    cb = ColorbarBase(ax_c,cmap=cmap,norm=norm,orientation='vertical',
                      label=cb_str)
    
    plt.savefig('Figures/'+title_str+'by State.png',bbox_inches='tight',dpi=200)
    plt.show()
    
if __name__ == "__main__":
    plt.ion()
    df = pandas.read_csv('openpv_all.csv',engine='python')

#    result = df['date_installed']
    inst_year = pandas.to_datetime(df['date_installed']).dt.year.values
    state = df['state'].str.slice(0,2).values
    size_kw = df['size_kw'].values
    cost = df['cost'].values
    cost_per_watt = df['cost_per_watt'].values
    pv_prod = df['annual_PV_prod'].values
    insol = df['annual_insolation'].values
    rebate = df['rebate'].str.replace(',','')
    rebate = rebate.str.replace('$','').values
    sales_tax = df['sales_tax_cost'].values
    tilt = df['tilt1'].values
    azimuth = df['azimuth1'].values
    reported_prod = df['reported_annual_energy_prod'].values
    pv_type = df['install_type'].str.lower().values
    states_count = df['state'].str.slice(0,2).value_counts()
    states_count = states_count.index.tolist()
    
    comm_index = numpy.where(pv_type == 'commercial')
    res_index = numpy.where(pv_type == 'residential')
    
#    type_str = ' Commercial'
#    type_index = comm_index
#    type_str = ' Residential'
#    type_index = res_index
    type_str = ' New'    
#    inst_year = inst_year[type_index]
#    state = state[type_index]
#    cost_per_watt = cost_per_watt[type_index]
#    cost = cost[type_index]
#    rebate = rebate[type_index]
#    size_kw = size_kw[type_index]
#    pv_prod = pv_prod[type_index]
#    reported_prod = reported_prod[type_index]
#    insol = insol[type_index]
    
    cpw_array = []
    net_cost = []
    net_cpw = []
    rebate = rebate.astype(numpy.float)
    prod_array = []
    insol_array = []
    rebate_array = []
    tot_rebate_array = []
    tot_cost_array = []
    tot_kw_array = []
    costperwatt_array = []
#    years = [1999.,2000.,2001.,2002.,2003.,2004.,2005.,2006.,2007.,2008.,2009.,2010.,2011.,2012.,2013.,2014.,2015.]
    years = [2015]
    for year in years:
        print year
        cpw_array = []
        net_cost = []
        net_cpw = []
        prod_array = []
        insol_array = []
        rebate_array = []
        tot_rebate_array = []
        tot_cost_array = []
        tot_kw_array = []
        costperwatt_array = []
        year_index = numpy.where(inst_year == year)
        year_state = state[year_index]
        year_cost_per_watt = cost_per_watt[year_index]
        year_cost = cost[year_index]
        year_rebate = rebate[year_index]
        year_size_kw = size_kw[year_index]
        year_pv_prod = pv_prod[year_index]
        year_reported_prod = reported_prod[year_index]
        year_insol = insol[year_index]
        for st in states_count:
            index = numpy.where(year_state == st)
            avg_cpw = numpy.nanmean(year_cost_per_watt[index])
            avg_cost = numpy.nanmean(year_cost[index])
            tot_cost = numpy.nansum(year_cost[index]/1000000000.)
            avg_rebate = max(0,numpy.nanmean(year_rebate[index]))
            avg_kw = numpy.nanmean(year_size_kw[index]*1000)
            tot_kw = numpy.nansum(year_size_kw[index]/1000.)
#            print avg_cost,avg_rebate
            tot_cost_array.append(tot_cost)
            tot_kw_array.append(tot_kw)
            costperwatt_array.append(tot_cost*1000./tot_kw)
            net_cost.append(avg_cost-avg_rebate)
            cpw_array.append(avg_cpw)
            net_cpw.append((avg_cost-avg_rebate)/avg_kw)
#            print avg_cpw, avg_cost/avg_kw
            avg_ann_prod = numpy.nanmean(year_pv_prod[index])
            avg_rep_prod = numpy.nanmean(year_reported_prod[index])
            prod_array.append(max(0,numpy.nanmean(year_reported_prod[index]/year_pv_prod[index])))
            avg_insol = numpy.nanmean(year_insol[index])
            insol_array.append(max(0,avg_insol))
            rebate_array.append(max(0,numpy.nanmean(100.*year_rebate[index]/year_cost[index])))
            tot_rebate_array.append(numpy.nansum(abs(year_rebate[index]))/1000000.)

#        state_names = StateConvert(states_count)
        state_names = states_count
    
#        print min(cpw_array),max(cpw_array)
        cpw_by_state = dict(zip(state_names,cpw_array))
#        print cpw_by_state
        net_cost_by_state = dict(zip(state_names,net_cost))
#        print min(net_cost),max(net_cost)
#        print net_cost_by_state

#        print 'tot_cost',min(tot_cost_array),max(tot_cost_array)
        tot_cost_by_state = dict(zip(state_names,tot_cost_array))
#        print tot_cost_by_state

#        print 'tot_kw',min(tot_kw_array),max(tot_kw_array)
        tot_kw_by_state = dict(zip(state_names,tot_kw_array))
#        print tot_kw_by_state
        
        costperwatt_by_state= dict(zip(state_names,costperwatt_array))
        print costperwatt_by_state
#        print'net_cpw',min(net_cpw),max(net_cpw)
        net_cpw_by_state = dict(zip(state_names,net_cpw))
#        print net_cpw_by_state
        
#        print'production %',min(prod_array),max(prod_array)
        production_by_state = dict(zip(state_names,prod_array))
#        print production_by_state
        
#        print'Insolation',min(insol_array),max(insol_array)
        insol_by_state = dict(zip(state_names,insol_array))
#        print insol_by_state
        
#        print'Rebate',min(rebate_array),max(rebate_array)
        rebate_by_state = dict(zip(state_names,rebate_array))
#        print rebate_by_state

        print'Rebate',min(tot_rebate_array),max(tot_rebate_array)       
        tot_rebate_by_state = dict(zip(state_names,tot_rebate_array))
        print tot_rebate_by_state

        BokehMap(tot_kw_by_state,0,1000,str(int(year))+type_str+' Photovoltaic Installations ', 'MW')
#        PlotMap(tot_kw_by_state,0,1000,str(int(year))+type_str+' Photovoltaic Installations ', 'MW')
#        PlotMap(tot_cost_by_state,0,5,str(int(year))+type_str+ ' Photovoltaic Investment ', '$B')
#        PlotMap(costperwatt_by_state,0,20,str(int(year))+type_str+' Photovoltaic Cost per Watt ', '$/W')
#        PlotMap(cpw_by_state,0,20,str(int(year))+' New Photovoltaic Cost per Watt ', '$/W')
        
#        PlotMap(insol_by_state,4,10,str(int(year))+' Average annual insolation ','Insolation')
#        PlotMap(rebate_by_state,0,100,str(int(year))+type_str+ ' Photovoltaic Rebate % ', '%')
#        PlotMap(tot_rebate_by_state,0,300,str(int(year))+type_str+ ' Photovoltaic Total Rebates ', '$M')
   
    