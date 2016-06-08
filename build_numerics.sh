#! /bin/bash

echo "BUILDING NUMERICS"

curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Info:\", \"color\": \"#00FF00\", \"text\":\"$2 completed succeffuly the openfpm_numerics test \" }] }" https://hooks.slack.com/services/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

