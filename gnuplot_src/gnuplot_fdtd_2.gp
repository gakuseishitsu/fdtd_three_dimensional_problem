#-------------------------------------------------------------------------------
# gnuplotの設定
#-------------------------------------------------------------------------------
reset
set nokey                 # 凡例の非表示
set xrange [0:0.6]       # x軸方向の範囲の設定
set yrange [0:0.6]         # y軸方向の範囲の設定
set zrange [-1:1]         # z軸方向の範囲の設定

set term gif animate delay 15      # 出力をgifアニメに設定
set output "fdtd_2.gif"  # 出力ファイル名の設定

#set contour
#set view 0,0
#unset surface
#set cntrparam levels 10

#set object 1 rect from 1.0, -0.4 to 1.1, 1 back linewidth 0 fillcolor rgb "#dcdcdc" fill solid 0.5
#set object 2 rect from 0, -0.4 to 0.008, 1 back linewidth 0 fillcolor rgb "#dcdcdc" fill solid 0.5
#set object 3 rect from 2, -0.4 to 2.016, 1 back linewidth 0 fillcolor rgb "#dcdcdc" fill solid 0.5

#-------------------------------------------------------------------------------
# 変数の設定
#-------------------------------------------------------------------------------
n0 = 1        # ループ変数の初期値
n1 = 70     # ループ変数の最大値
dn = 1  # ループ変数の増加間隔

#-------------------------------------------------------------------------------
# ループの開始
#-------------------------------------------------------------------------------
load "gnuplot_fdtd_2.plt" 