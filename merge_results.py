import pandas as pd
import os
import re

def merge_deseq_results():
    # 设置结果目录和输出文件路径
    result_dir = 'DESeq2_results'
    output_file = 'merged_deseq_results.csv'
    
    # 获取所有CSV文件
    csv_files = [f for f in os.listdir(result_dir) if f.endswith('.csv')]
    if not csv_files:
        print(f'错误: 在{result_dir}目录中未找到CSV文件')
        return
    
    # 初始化合并结果DataFrame
    merged_df = None
    
    for file in csv_files:
        # 从文件名提取KO基因名称 (例如从'BCL11A_vs_Nontargeting_control.csv'提取'BCL11A')
        ko_gene = re.sub(r'_vs_Nontargeting_control\.csv$', '', file)
        file_path = os.path.join(result_dir, file)
        
        try:
            # 读取CSV文件，只选择需要的列
            df = pd.read_csv(
                file_path,
                usecols=['gene_id', 'log2FoldChange', 'padj'],
                dtype={'gene_id': str}
            )
        except Exception as e:
            print(f'读取文件{file}时出错: {str(e)}')
            continue
        
        # 重命名列，添加KO基因前缀
        df.rename(columns={
            'log2FoldChange': f'{ko_gene}_log2FoldChange',
            'padj': f'{ko_gene}_padj'
        }, inplace=True)
        
        # 合并到主DataFrame
        if merged_df is None:
            merged_df = df
        else:
            merged_df = merged_df.merge(df, on='gene_id', how='outer')
    
    # 保存合并结果
    if merged_df is not None:
        merged_df.to_csv(output_file, index=False)
        print(f'成功合并{len(csv_files)}个结果文件')
        print(f'合并结果已保存至: {output_file}')
        print(f'表格维度: {merged_df.shape[0]}个基因 × {merged_df.shape[1]-1}个指标')
    else:
        print('未合并任何数据')

if __name__ == '__main__':
    merge_deseq_results()

