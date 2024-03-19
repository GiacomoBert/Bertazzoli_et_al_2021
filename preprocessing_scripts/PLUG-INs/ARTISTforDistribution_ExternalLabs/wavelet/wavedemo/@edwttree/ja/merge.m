%MERGE   各ノードのデータをマージします(組み替えます)
%   X = MERGE(T,N,TNDATA) は、ノード N の子ノードで関連づけられたデータ
%   を使って、データツリー T のノード N で関連づけられたデータ X を組み
%   替えます。
%
%   TNDATA は、TNDATA{k}が N の k 番目の子ノードとなるよう関連づけられた%   データを含むような、(ORDER x 1)行(1 x ORDER)列となるセル配列です。
%
%   IDWT を用います。
 
%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Sep-1999.


%   Copyright 1995-2003 The MathWorks, Inc.
