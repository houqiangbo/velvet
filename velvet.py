#pandas 对 fasta 文件进行分组
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from functools import reduce
import pandas as pd

#试试直接分组，用parse,把 id/seq分两个列表

dict1 = dict()
#name = []
set1 = set()
for ret in SeqIO.parse(open("/home/hou/htf/getall.fasta","r"),"fasta"):
    if 'R' in str(ret.id[-3:]) or 'r' in str(ret.id[-3:]):
        dict1.update({ret.id : str(ret.seq.reverse_complement())})
    else:
        dict1.update({ret.id : str(ret.seq)})
    set1.add(ret.id.split("_")[0] + "_" + ret.id.split("_")[1])
     #    name.append(ret.id.split("_")[0] + "_" + ret.id.split("_")[1])
#print(name)

#print(len(set1))
dict2  = dict()
for i in set1:
    seq = ""
    for j in dict1:
        if i.split("_")[0] + "_" + i.split("_")[1] == j.split("_")[0] + "_" + j.split("_")[1]:
            seq += dict1[j] + "|"
    dict2.update({i.split("_")[0] + "_" + i.split("_")[1] : seq})

#readlist=[]
#for key in dict2:
#        readlist=dict2[key].split("|")[0:-1]

def get_weight(s1,s2):               #通过两条序列的overlap计算出权值
    l = min(len(s1),len(s2))
    while l>0:
        if s2[:l] == s1[-l:]:
            return l
        else:
            l-=1
    return 0

def print_result(s1,s2):            #将两条序列去除首尾overlap后合并
    weight = get_weight(s1,s2)
    s = s1 + s2[weight:]
    #print(s)
    return s

def dir_graph(l,t=3):              #得到满足条件的有向图(权值大于等于t)
    graph = {}
    for i in l:
        VW = []
        for j in l:
            if i!=j:
                weight = get_weight(i,j)
                if weight >= t:
                    VW.append(j)
        graph[i] = VW
    #print(graph)
    for i in graph.keys():        #不能有孤立顶点
        if not graph[i]:
            break
            count = get_in_V(graph,i)
            if count ==0:
                graph.clear()
                print('The sequence:\n"{0}"\n can\'t align with others!'.format(i))
                break
    return graph

def get_in_V(graph,v):                   #得到某顶点入度
    count = 0
    all_in = reduce(lambda x,y:x+y,graph.values())
    for i in all_in:
        if i == v:
            count+=1
    return count

def aligner(graph,topo=[]):             #得出顶点顺序
    while graph:
        V = graph.keys()
        for i in V:
            flag = 1
            in_num = get_in_V(graph,i)
            if in_num ==0:
                topo.append(i)
                graph.pop(i)
                flag = 0
                break
        if flag:                        #存在环
            #print('The t score is too small!')
            return None
        else:
            aligner(graph,topo)
    return topo



li = []
velvet_dict = {}
for key in dict2:
#    l = []
    readlist = dict2[key].split("|")[:-1]
#    print(readlist[:-1])
    print("\n")
    print(key)
#    pass
#print(dict2)
       
    graph = dir_graph(readlist,t=3)
#    l.append(graph)
    topo = aligner(graph)
    if topo:
        result = reduce(print_result,topo)
    else:
        result = topo
#    l.append(result)
#    li.append(str(graph) + "|" + str(result))
#    print(key)
    print(result)
    velvet_dict.update({key : result})
#print(li)
#print(velvet_dict)
#print(len(dict2))

#for key in velvet_dict:
#    print(key)

#    print(result)
#        print(">"+key+"\n"+result)
#    except TypeError:
#        continue
#        print(key + "error")
#    else:
#        print(">"+key+"\n"+result)
        

#print(readlist)


