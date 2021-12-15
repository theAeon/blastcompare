#!/usr/bin/env pwsh
$file=$args[0]
$dir=Get-ChildItem $file -File -Name
$dir | ForEach-Object -ThrottleLimit 8 -parallel {
    $a=Import-Csv -Delimiter "`t" -Path $using:file/$_ -Header name, id
    $a | Select-Object -Property id | Format-Table -HideTableHeaders | Out-File -FilePath ./fastas/$_.ids
    blastdbcmd -db combined -entry_batch ./fastas/$_.ids -out ./fastas/$_.fasta 2>> ./fastas/$_.err
}
