% MERGE   各ノードのデータをマージします(組み替えます)
%
%   X = MERGE(T,N,TNDATA) は、ノード N の子ノードで関連づけられたデータ
%   を使って、データツリー T のノード N で関連づけられたデータ X を組み
%   替えます。
%
%   TNDATA は、TNDATA{k}が N の k 番目の子ノードとなるよう関連づけられた
%   データを含むような、(ORDER x 1) または (1 x ORDER) となるセル配列です。
%
%   注意: 
%   このメソッドは、オブジェクトのクラスに対して多重定義されています。
%   ----------
%   クラス DTREE に対して、MERGE は、N の最も左側の子ノードのデータを
%   割り当てます。
%   ----------


%   Copyright 1995-2004 The MathWorks, Inc.
