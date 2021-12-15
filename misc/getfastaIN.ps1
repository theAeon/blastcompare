#!/usr/bin/env pwsh
$file=$args[0]
$dir=Get-ChildItem $file -File -Name
$dir | ForEach-Object -ThrottleLimit 8 -parallel {
    Import-Csv -Delimiter "`t" -Path $using:file/$_ -Header name, id | Group-Object -Property name | Format-Table -HideTableHeaders | Out-File -FilePath ./inseq/$_.ids
    blastdbcmd -db combined -entry_batch ./inseq/$_.ids -out ./inseq/$_.fasta 2>> ./inseq/$_.err

}
