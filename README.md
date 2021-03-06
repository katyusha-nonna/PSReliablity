## 发电系统可靠性评估
### author: Katyusha

用于计算发电系统或互联发电系统的发电裕度表&可靠性指标。

使用程序前，使用者需按照输入数据的模板组织输入数据，格式如下：

**格式：**

[模式选择]	 [分析步长] [分析周期] [元件总数] [发电机数] [负荷数] [输电线路数] [变压器数] 

[发电机名称] [发电机ID] [发电机容量] [发电机故障率lambda] [发电机修复率miu]

[负荷名称] [负荷ID] [负荷容量]
[负荷时序曲线]

[线路名称] [线路ID] [线路容量] [线路故障率lambda] [线路修复率miu]

[变压器名称] [变压器ID] [变压器容量] [变压器故障率lambda] [变压器修复率miu]

[系统拓扑/计算表达式]

**说明：**
* [1]	可针对包含四种元件的发电系统或互联系统的可靠性进行评估。其中各个元件的条目数量要和第一行的基本数目参数对应。且各元件的输入顺序应该保持；发电机-负荷-线路-变压器。总元件数目包含了发电机、负荷、输电线和变压器，值应该与上述四项参数之和相等。
* [2]	负荷的时序曲线离散点数目应该和第一行的控制参数中分析周期一致。发电机、线路与变压器的故障率和修复率需统一单位至 次/日。
* [3]	系统的拓扑的表示采用中缀计算表达式，直接根据中缀表达式的形式给出计算流程即可。串联等效计算请用运算符*，并联卷积计算请用运算符+，优先级顺序为：括号("()") > 串联("*") > 并联("+")。
* [4]	由于采用了灵活输入，模式选择暂时不起作用，只作为预留的接口。
* [5]	输入文件需命名为input.txt。

**举例：**

输入文件“input.txt”和输出结果"output.txt"已放置在Example文件夹内

系统拓扑如图所示

![](./Image/system.png)




