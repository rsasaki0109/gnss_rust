# 次に比較するSD ambiguityの共通オフセット

固定した元C++の[src/algorithms/rtk_bias.cpp](https://github.com/rsasaki0109/gnssplusplus-library/blob/72f3b7c2c4c088dc575286a6f0edf4e407bc499a/src/algorithms/rtk_bias.cpp#L571)
では、衛星系ごとに現在のb=L_SD−P_SD/λと既存のSD ambiguity状態との差を平均し、
保持している同じ系のSD状態にその平均を加えてから、新しい状態をbで作る。
この状態の単位はcycleで、初期分散はrtk_state.cppで900 cycle²。単位変換の差ではない。

全ての保持状態への共通平行移動は、その状態同士のDDを変えない。
一方、新しく加わる状態の平均と保持状態のgaugeをそろえる点で意味がある。
現在のRust priorは保持状態をそのまま残し、新しい状態を直接bで初期化している。
この差を、欠測/再捕捉・基準変更・コードの共通偏り・ゼロの有効状態を含む独立ケースで
検証する。PPC外れの原因と確定したわけではなく、この記録の実行へ遡って適用していない。

既存の数値許容差を変えず、適用する場合は新しい実行記録として比較する。
