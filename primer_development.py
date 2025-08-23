import streamlit as st
import math
import itertools

st.title("FASTA対応診断用プライマー候補配列探索ツール")

iupac_dict = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'C', 'G', 'T']): 'N'
}

iupac_reversed_dict = {v: k for k, v in iupac_dict.items()}

def get_iupac(bases):
    return iupac_dict.get(frozenset(bases), 'N')

def calc_base_frequencies(sequences, pos):
    bases = [seq[pos] for seq in sequences if seq[pos] in ['A', 'C', 'G', 'T']]
    total = len(bases)
    if total == 0:
        return {}
    freq = {}
    for b in ['A', 'C', 'G', 'T']:
        freq[b] = bases.count(b) / total
    return freq

def calc_gc_content(seq):
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    length = len(seq.replace('-', ''))
    if length == 0:
        return 0
    return (gc_count / length) * 100
    
def calculate_tm_with_iupac(seq):
    
    parameter = {
    'AA':[-7.9, -22.2],
    'TT':[-7.9, -22.2],
    'AT':[-7.2, -20.4],
    'TA':[-7.2, -21.3],
    'CA':[-8.5, -22.7],
    'TG':[-8.5, -22.7],
    'GT':[-8.4, -22.4],
    'AC':[-8.4, -22.4],
    'CT':[-7.8, -21.0],
    'AG':[-7.8, -21.0],
    'GA':[-8.2, -22.2],
    'TC':[-8.2, -22.2],
    'CG':[-10.6, -27.2],
    'GC':[-9.8, -24.4],
    'GG':[-8.0, -19.9],
    'CC':[-8.0, -19.9]
        }
    
    expanded = [iupac_reversed_dict[base] for base in seq]
    expanded_seq = [''.join(codon) for codon in itertools.product(*expanded)]
    tm_range = []

    for i in expanded_seq:
        pair = [i[j:j+2] for j in range(len(i)-1)]
        sumH = sum([parameter[i][0] for i in pair]) 
        sumS = sum([parameter[i][1] for i in pair])
        gc_content = calc_gc_content(i)/100
        if i[0] in ['G', 'C']:
            ini_H = 0.1
            ini_S = -2.8
        elif i[0] in ['A', 'T']:
            ini_H = 2.3
            ini_S = 4.1
        tm_format = (1000*(sumH+ini_H)/(ini_S+sumS+1.987*math.log(primer_molar)))
        ln_m = math.log(salt_molar)
        correction = ((4.29 * gc_content - 3.95) * ln_m + 0.94 * (ln_m ** 2)) * 1e-5
        tm_corrected = round(tm_format/(1 + correction*tm_format) - 273.15, 1)
        tm_range.append(tm_corrected)

    return [min(tm_range), max(tm_range)]


def get_base_frequencies_at_positions(sequences, start, end):
    position_frequencies = []
    for i in range(start, end):
        freq = calc_base_frequencies(sequences, i)
        position_frequencies.append(freq)
    return position_frequencies

def read_fasta_alignment(lines):
    sequences = []
    seq_names = []
    current_seq = []

    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_seq:
                sequences.append(list(''.join(current_seq)))
                current_seq = []
            seq_names.append(line)
        else:
            current_seq.append(line)
    if current_seq:
        sequences.append(list(''.join(current_seq)))
    return seq_names, sequences

# ユーザーが閾値をスライダーで設定可能に
st.sidebar.header("プライマー条件設定")
min_tm = st.sidebar.slider("最小Tm値 (℃)", 0, 100, 55)
max_tm = st.sidebar.slider("最大Tm値 (℃)", 0, 100, 65)
min_gc = st.sidebar.slider("最小GC含有率 (%)", 0, 100, 40)
max_gc = st.sidebar.slider("最大GC含有率 (%)", 0, 100, 60)
min_len = st.sidebar.slider("最小塩基長", 0, 100, 20)
max_len = st.sidebar.slider("最大塩基長", 0, 100, 30)
fr = st.sidebar.slider("完全一致率", 0, 100, 90)
primer_molar = st.number_input('primer濃度(nM)', value=500) * 1e-9
salt_molar = st.number_input('塩濃度(mM)', value=50) * 1e-3

def analyze_block(sequences, block_num=1):
    if len(sequences) == 0:
        st.write(f"--- 配列なし ---")
        return

    length_set = set(len(seq) for seq in sequences)
    if len(length_set) != 1:
        st.warning(f"配列長が揃っていません。処理を中止します。")
        return

    seq_len = length_set.pop()
    st.write(f"### 配列数={len(sequences)}, 配列長={seq_len}")

    max_rates = []
    consensus_seq = []

    for i in range(seq_len):
        freq = calc_base_frequencies(sequences, i)
        if not freq:
            consensus_seq.append('N')
            max_rates.append(0)
            continue
        max_freq = max(freq.values())
        max_rates.append(max_freq)
        if max_freq == 1.0:
            # 完全一致は単一塩基
            base = max(freq, key=freq.get)
            consensus_seq.append(base)
        else:
            # 完全一致でない場合は存在する塩基すべてで混合塩基
            bases_present = [b for b, r in freq.items() if r > 0]
            consensus_seq.append(get_iupac(bases_present))

    consensus_str = ''.join(consensus_seq)

    threshold = 1.0

    candidates = []
    progress_bar = st.progress(0)

    for window_size in range(min_len, max_len + 1):
        for start in range(seq_len - window_size + 1):
            window_rates = max_rates[start:start + window_size]
            high_rate_count = sum(r >= threshold for r in window_rates)
            if high_rate_count >= window_size * fr * 1e-2:
                # ギャップ含む配列除外
                has_gap = any('-' in seq[start:start + window_size] for seq in sequences)
                if has_gap:
                    continue

                primer_seq = consensus_str[start:start + window_size]
                tm = calculate_tm_with_iupac(primer_seq)
                gc = calc_gc_content(primer_seq)
                fullmatch_count = sum(1 for rate in window_rates if rate == 1.0)*100/len(primer_seq)

                if min_tm <= tm[0] and tm[1] <= max_tm and min_gc <= gc <= max_gc and fr <= fullmatch_count:
                    candidates.append((start + 1, start + window_size, primer_seq, tm, gc, fullmatch_count))
        progress_bar.progress(min(((window_size - min_len + 1)/(max_len - min_len + 1), 1.0)))

    if candidates:
        st.subheader(f"プライマー候補領域（開始-終了 : 配列 (Tm℃, GC%, 完全一致率)）")
        sorted_candidates = sorted(candidates, key=lambda x: x[0])
        for start_pos, end_pos, seq, tm, gc, fullmatch_count in sorted_candidates:
            st.text(f"{start_pos}-{end_pos}: {seq} (Tm={tm[0]}-{tm[1]}℃, GC={gc:.1f}%, 完全一致率={fullmatch_count:.1f}%)")
            # 混合塩基がある位置の塩基割合を表示
            freqs = get_base_frequencies_at_positions(sequences, start_pos - 1, end_pos)
            iupac_positions = []
            for i, base in enumerate(seq):
                if base not in ['A', 'T', 'G', 'C']:
                    freq = freqs[i]
                    if freq:
                        detail = ', '.join([f"{b}...{round(r*100,1)}%" for b, r in freq.items()])
                        iupac_positions.append(f"{base} : {detail}")
            if iupac_positions:
                for detail in iupac_positions:
                    st.write(f"　　  {detail}")
            
    else:
        st.write(f"条件に合うプライマー候補領域が見つかりませんでした。")

uploaded_file = st.file_uploader("FASTA形式アライメントファイルをアップロード（utf-8）", type=["fa", "fasta", "txt"])

if uploaded_file is not None:
    lines = uploaded_file.read().decode('utf-8').splitlines()
    seq_names, sequences = read_fasta_alignment(lines)

    if len(sequences) == 0:
        st.error("配列が読み込めませんでした。FASTA形式ファイルを確認してください。")
    else:
        length_set = set(len(seq) for seq in sequences)
        if len(length_set) != 1:
            st.error("配列長が揃っていません。アラインメントされたFASTAファイルをアップロードしてください。")
        else:
            with st.status("処理中...", expanded=True) as status:
    # 処理中の内容など表示したければここに書く
                analyze_block(sequences, 1)
                status.update(label="完了しました", state="complete")
