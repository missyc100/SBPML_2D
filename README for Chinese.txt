------------------------------------------------------------快速地运行我们的案例---------------------------------------------------------------------------
1、修改U1004-E24.for中的读取文件路径（参阅U1004-E24.for文件第13~17行）
2、在ABAQYS中建立Job-PML.inp文件的job并关联U1004-E24.for
-----------------------------教程分为两个部分：第1部分是SBPML UEL单元，第2部分是域缩减法----------------------------------------------------------
我使用的环境是ABAQUS 2022, MATLAB 2022a, Intel(R) MPI Library 2019 Update 9 
Calculated wave velocity文件夹的主程序是main.m
Generate matching layer parameters with one click文件夹的主程序是main1.m
DRM文件夹的主程序是main.m

1、在MATLAB中安装Abaqus2Matlab.mlappinstall插件
2、运行Calculated wave velocity\main得到建议的单元尺寸、匹配层厚度、吸收系数
3、 建立类似model1.cae的模型，SBPML域需要在PART1中创建，SBPML域每一层都必须保证平行，分析步选择隐式动力，矩阵储存方式选择非对称
4、 在PART1中建立集合
cross			SBPML单元和CPE4单元交界面处的节点集
element_all		PART1所有单元集合
element_set1		CPE4单元集合
element_set2		SBPML单元集合
fixnode			固定节点集
node_all		PART1所有节点集合
node_set1		CPE4节点集合
node_set2		SBPML节点集合
PML1			材料1的SBPML单元集合
PML2			材料2的SBPML单元集合
m1			材料1的CPE4单元集合
m2			材料2的CPE4单元集合
5、在cae中输出INP文件，将对应名称的集合复制到	'\Generate matching layer parameters with one click\input'	路径中
6、修改\Generate matching layer parameters with one click\main1.m中120行匹配层层数,152行材料数目, 教程中是2种材料
7、运行main1.m
8、在cae文件所在目录创建read文件夹，将第7步中输出的output1/element_set1.txt,output1/element_set2.txt和output3中的所有文件复制到read文件夹
9、修改第5步输出的INP文件
	具体步骤是	(1) 删除PART1中的单元连通性
			(2) 将UEL定义的格式复制到相应位置（参阅Job-PML.inp   813行）
			(3) 修改CPE4赋予的材料属性   (参阅Job-PML.inp   1092、1095行)
10、修改U1004-E24.for中的读取文件路径，单元节点个数，材料属性和吸收参数（参阅U1004-E24.for文件第2行）
11、在ABAQYS中建立Job-PML.inp文件的job并关联U1004-E24.for
-----------------------------------------至此SBPML可以正常工作，以下是地震动输入部分---------------------------------------------------------------------
1、前述第5步中输出的INP文件，将对应名称的集合复制到	'\DRM\input'	路径中
2、选择一个PART进行DRM输入，这个例子选择的是PART1，将PART1输出到一个新的INP文件，并在INP文件最后加入输出质量、刚度矩阵的命令（参阅Job-KKMM.inp   1641~1650行）
3、在ABAQUS中运行第2步得到的INP文件，运行结束会得到两个*.mtx文件，将这两个文件复制到   '\DRM\MTX'   路径中
4、确定想要输入的地震动时程（参阅   '\DRM\main.m'   第63行）
5、确定想要输入的波的类型，可选1、P波，2、SV波（参阅   '\DRM\main.m'   第72行）。地震动入射角度（参阅   '\DRM\main.m'   第65行）。注意Fei1可调节，Fei2为0.
6、修改基岩材料属性（参阅   '\DRM\OneDimension_FE_AB.m'   第8~10行）
7、修改地层材料属性（参阅   '\DRM\material_constant_Solid.txt'   ）
	注意：	需要查看'\DRM\main.m'运行到72行时矩阵YY的内容，这会告诉你需要定义多少层材料
		不要删除'\DRM\material_constant_Solid.txt' 第一行的lamd内容
8、在cae文件所在目录创建case1文件夹，将\DRM\output1,\DRM\output2中的内容以及\DRM\output4文件夹复制到case1文件夹（参阅\DRM\case1）
9、将第8步中得到的txt文件include到INP文件中（参阅Job-PML-DRM.inp文件，第2404、2446、2476行）
10、在ABAQYS中建立Job-PML-DRM.inp文件的job并关联U1004-E24.for
---------------------------------------------------------结束-----------------------------------------------------------------------------------------------------------